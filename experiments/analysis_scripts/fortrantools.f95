
subroutine my_lindemann(ref_coords, coords, num_nb, list_nb,&
     box, cutoff, MAXnb, N, D, lind)
  !===========================================================
  ! UNFINISHED I THINK! Do not use this as is in the code.
  !===========================================================
  integer, intent(in):: N, D, MAXnb
  real(8), intent(in)   :: coords(N,D), ref_coords(N,D)
  real(8)               :: delta(N,N)
  real(8), intent(in)   :: box(D), cutoff
  integer, intent(in):: num_nb(N), list_nb(N, MAXnb)
  real(8), intent(out)  :: lind
  integer            :: i, j, k
  real(8)              :: rsq, dis, r(D)
  integer            :: cnt
  
  lind = 0.
  delta = coords-ref_coords
  
  do i=1, N-1
     do k=1, num_nb(i)
        j = list_nb(i, k)
        r = delta(i,:) - delta(j,:)
        call PBC(r , box)
        lind = lind + sum(r*r)
        cnt = cnt + 1
     end do
  end do
  
  if (cnt>0) then
     lind = lind / cnt
  else
     lind = 0.
  endif
  print*, 'all done!', lind
  return 
end subroutine my_lindemann


subroutine distance_matrix(coords, psi6_re, psi6_im, mus,&
     local_rho, box, cutoff, gsq_flag, rho_0,MAXnb, N, D,& 
     nbins, nbinsq, distance, distance_x,distance_y,&
     num_nb,list_nb, g, g6, g6_re, g6_im, sq, g_ori,&
     g_dp, g_dp_tr, g_pr, s_pr)
  !===========================================================
  ! This function computes lots of things at the same time
  ! for efficiency.
  !
  ! Compute distance(N,N), which is the pairwise distance
  !  between all of the atoms.
  !
  ! Compute distance_x(N,N), distance_y(N,N) which are just
  !  x and y components of distance between pairs.
  !
  ! Compute list_nb(N,MAXnb) for which each atom tracks
  !  a list of its neighbours.
  !
  ! Compute a num_nb(N), a count for the number of neighbours
  !  each atom has.
  !
  ! Compute g(nbins,5), which is the pairwise correlation
  !  function, and the second index (:,i) with i = 1,...,5
  !  correspond to different correlation functions:
  !  i = 1 is just the distances (r in g(r)),
  !  i = 2 is the full pairwise correlation function,
  !  i = 3 is the pairwise correlation function between only
  !  neighbours in a dilute region,
  !  i = 4 is the pairwise correlation function between
  !  only neighbours in a dense region, and
  !  i = 5 is the pairwise correlation function where one
  !  atom is in the dense region and the other is in the
  !  dilute region.
  !
  ! Compute g6(nbins,5), which is modulus of the hexagonal
  !  order parameter correlation function, <psi_i^{*}psi_j>.
  !
  ! Compute g6re(nbins,5) and g6im(nbins,5), which are the
  !  real and imaginary parts of the hexagonal order
  !  parameter correlation function.
  !
  ! Compute sq(nbinsq,nbinsq,3) which is the structure
  !  factor of the system. the third index (:,:,i) with
  !  i = 1,2,3 correspond to:
  !  i = 1 is the wavevector q's x component, q_x,
  !  i = 2 is the wavevector q's y component, q_y, and
  !  i = 3 is the structure factor.
  !
  ! Compute g_ori(nbins,5), which is the orientational
  !  correlation function < u(r_i)u(r_j) >_{neighbours}.
  !
  ! Compute g_dp(nbins,5), which is
  !  <\hat{u}_{ij}*\hat{r}_{ij}>_{neighbours} at each
  !  position.
  !
  ! Compute g_dp_tr(nbins,5), which is
  !  <\hat{u}_{ij}
  !  -(\hat{u}_{ij}*\hat{r}_{ij})\hat{r}_{ij}>_{neighbours}.
  !
  ! Compute g_pr(nbins,5), which is
  !  <u_{ij}*r_{ij}>_{neighbours}.
  !
  ! Compute s_pr(nbins,5), which is
  !  <u_{ij}*r_{ij}>_{neighbours,atoms}.
  !
  ! N is number of atoms in the simulation.
  ! D is dimension.
  ! MAXnb is maximum number of neighbours expected (100 is
  !  a safe bet).
  ! nbins is the number of bins to divide r in correlation
  !  function up into.
  ! nbinssq is the number of bins to divide each axis of the
  !  (2D) q vector for the structure factor into.
  ! gsq_flag = 1 if you want to compute structure factor and
  !  correlation functions along with distance matrices.
  ! coords(N,D) is the coordinates of all atoms.
  ! psi_6re(N), psi6_im(N) are the real and imaginary parts
  !  of the hexagonal order parameter for each atom.
  ! mus(N,D) are the dipole vectors for each atom.
  ! local_rho(N) is the local density around each atom.
  ! box(D) is the dimensions of the periodic box.
  ! cutoff is the distance at which two atoms are no
  !  longer neighbours.
  ! rho_0 is the average density of the system.
  !===========================================================
  implicit none
  integer, intent(in)  :: N, D, MAXnb, nbins, nbinsq, gsq_flag
  real(8), intent(in)  :: coords(N,D), psi6_re(N), psi6_im(N)
  real(8), intent(in)  :: mus(N,D), local_rho(N)
  real(8), intent(in)  :: box(D), cutoff, rho_0
  real(8), intent(out) :: distance(N, N)
  real(8), intent(out) :: distance_x(N, N),distance_y(N, N)
  integer, intent(out) :: list_nb(N, MAXnb), num_nb(N)
  real(8), intent(out) :: g(nbins, 5), g6(nbins,5), g6_re(nbins,5)
  real(8), intent(out) :: g6_im(nbins,5), sq(nbinsq,nbinsq,3)
  real(8), intent(out) :: g_ori(nbins,5), g_dp(nbins,5)
  real(8), intent(out) :: g_dp_tr(nbins,5), g_pr(nbins,5), s_pr(2:5)
  integer              :: i, j, ix, k, num_q_vec, kx, ky
  integer              :: ij_type, Ndummy(2:5)
  real(8)              :: rsq, dis, r(D), cutoff_sq, dr
  real(8)              :: pi=3.14159265359, dp, pij_hat(D), dummy
  real(8)              :: re, im, dq, q(D), theta, qmin, qmax, rq
  real(8)              :: sin_rq(nbinsq,nbinsq), cos_rq(nbinsq,nbinsq)
  

  CALL RANDOM_SEED
  num_q_vec = 1
  distance = 0.0
  num_nb = 0
  cutoff_sq = cutoff *  cutoff 
  
  dr =  dsqrt(box(1)*box(2)/4.)  / (nbins - 1)
  
  qmin = -20 !2*pi/dsqrt(box(1)*box(2)/4.) / 2
  qmax = 20.
  dq =  (qmax-qmin) / (nbinsq - 1)
  g = 0
  g6_re = 0
  g6_im = 0
  g6 = 0
  sq = 0
  g_ori = 0
  g_dp = 0
  g_dp_tr = 0
  g_pr = 0
  s_pr = 0
  sin_rq = 0
  cos_rq = 0
  Ndummy=0


  do i=1, N-1
     ! sum over atoms 1 to N-1     
     if (gsq_flag == 1) then
        ! sq
        r = coords(i,:)
        call PBC(r,box)
        ! compute sum_{r} e^{i q * r} for each q point
        do kx=1, nbinsq
           do ky=1, nbinsq
              q(1) = qmin + (kx+0.5)*dq
              q(2) = qmin + (ky+0.5)*dq
              rq = dot_product(r, q)
              sin_rq(kx,ky) = sin_rq(kx,ky) + dsin(rq)
              cos_rq(kx,ky) = cos_rq(kx,ky) + dcos(rq)
           enddo
        enddo
     endif


     do j=i+1, N
        ! combined with above sum over 
        ! sum over each pair of atoms
        r = coords(i,:) - coords(j,:) 
        call PBC(r,box) 
        rsq = sum(r*r)
        dis = dsqrt(rsq)
        if (rsq <= cutoff_sq) then
           ! increase neighbour count for both i and j atoms
           num_nb(i) = num_nb(i) + 1
           num_nb(j) = num_nb(j) + 1
           ! track where neighbours are in list
           list_nb(i, num_nb(i)) = j
           list_nb(j, num_nb(j)) = i
        end if
        ! track neighbour matrix for distance
        distance(i,j) = dis
        distance(j,i) = dis
        ! track neighbour matrix for positions
        distance_x(i,j) = r(1)
        distance_x(j,i) = -r(1)      
        distance_y(i,j) = r(2)
        distance_y(j,i) = -r(2)
        
        
        if ((local_rho(i)<rho_0) .AND. (local_rho(j)<rho_0)) then
           ! in dilute region
           ij_type = 3
        else if ((local_rho(i)>=rho_0) .AND. (local_rho(j)>=rho_0)) then
           ! in dense region
           ij_type = 4
        else
           ! in boundary region
           ij_type = 5
        end if
        
        
        ! count total pairs
        Ndummy(2)=Ndummy(2)+1
        ! count number of pairs in dense, dilute, or boundary
        Ndummy(ij_type)=Ndummy(ij_type)+1
        
        
        

        if (gsq_flag == 1) then 
           ! g_of_r
           ! integer number of dr spacings
           ix = int(dis/dr)
           if (ix <= nbins) then
              ! g counts the number of atoms of each type, at a specified
              ! bin number (i.e. distance from each other), dN
              ! it will eventually be the pair correlation function g(r)
              ! defined as <rho> g(r) dr^2 = dN/(N-1).
              g(ix, 2) = g(ix, 2) + 2 
              g(ix, ij_type) = g(ix, ij_type) + 2
              re =   ( psi6_re(i) * psi6_re(j)  +  psi6_im(i) * psi6_im(j) )
              im =   ( -psi6_re(i) * psi6_im(j) +  psi6_im(i) * psi6_re(j) )
              ! g6 is hexagonal ordering pair correlation, defined as
              !  <psi(r_i)*psi(r_j)>_{neighbours} where the sum is over
              ! all neighbours within cutoff distance
              g6_re(ix, 2) = g6_re(ix, 2) + 2 * re    
              g6_re(ix, ij_type) = g6_re(ix, ij_type) + 2 * re    
              g6_im(ix, 2) = g6_im(ix, 2) + 2 * im   
              g6_im(ix, ij_type) = g6_im(ix, ij_type) + 2 * im   
              
              dummy = 2*sum(mus(i,:)*mus(j,:))
              ! g_ori is orientational ordering at each distance
              ! < u(r_i)u(r_j) >_{neighbours}
              g_ori(ix, 2) = g_ori(ix, 2) + dummy
              g_ori(ix, ij_type) = g_ori(ix, ij_type) + dummy
              
              pij_hat = (mus(i,:)-mus(j,:))
              dummy = 2*sum(pij_hat*r)
              ! g_pr is <u_{ij}*r_{ij}>_{neighbours}
              g_pr(ix, 2) = g_pr(ix, 2) + dummy
              g_pr(ix, ij_type) = g_pr(ix, ij_type) + dummy

              ! s_pr is u_{ij}*r_{ij} averaged over all atoms, i.e.
              ! <u_{ij}*r_{ij}>_{neighbours,atoms}
              s_pr(2) = s_pr(2) + sum(pij_hat*r)
              s_pr(ij_type) = s_pr(ij_type) + sum(pij_hat*r)
              
              ! g_dp is \hat{u}_{ij}*\hat{r}_{ij} at each position
              pij_hat = pij_hat / dsqrt(sum(pij_hat * pij_hat))
              dp = sum(pij_hat*r)/dis 
              g_dp(ix, 2) = g_dp(ix, 2) +  2*dp  
              g_dp(ix, ij_type) = g_dp(ix, ij_type) +  2*dp  

              ! g_dp_tr is
              ! \hat{u}_{ij}-(\hat{u}_{ij}*\hat{r}_{ij})\hat{r}_{ij}
              pij_hat = pij_hat - dp*r/dis
              dummy = 2*dsqrt(sum(pij_hat*pij_hat))
              g_dp_tr(ix, 2) = g_dp_tr(ix, 2) + dummy
              g_dp_tr(ix, ij_type) = g_dp_tr(ix, ij_type) + dummy 
              
           end if

        end if
        
     enddo
  enddo
  
  
  if (gsq_flag .eq. 1) then 
     
     
     do i=1, nbins
        ! store location of bin centre   
        g(i,1) = (i+0.5)*dr
        g6(i,1) = (i+0.5)*dr
        g6_re(i,1) = (i+0.5)*dr
        g6_im(i,1) = (i+0.5)*dr
        g_ori(i,1) = (i+0.5)*dr
        g_dp(i,1) = (i+0.5)*dr
        g_dp_tr(i,1) = (i+0.5)*dr
        g_pr(i,1) = (i+0.5)*dr
        
        do ij_type = 2, 5
           if (g(i, ij_type)>0) then
              ! normalize things by particle number in the shell, dN
              ! where dN is currently g
              g6_re(i,ij_type) = g6_re(i,ij_type) / g(i, ij_type) 
              g6_im(i,ij_type) = g6_im(i,ij_type) / g(i, ij_type)
              g6   (i,ij_type) = dsqrt (g6_re(i,ij_type)*g6_re(i,ij_type)&
                   + g6_im(i,ij_type)*g6_im(i,ij_type) )
              g_ori(i,ij_type) = g_ori(i,ij_type) / g(i, ij_type) 
              g_dp(i,ij_type) = g_dp(i,ij_type) / g(i,ij_type) 
              g_dp_tr(i,ij_type) = g_dp_tr(i,ij_type) / g(i,ij_type)
              g_pr(i,ij_type) = g_pr(i,ij_type) / g(i,ij_type)
              ! since N_pair = N(N-1)/2, where N is number of particles
              ! in the system, and since N_pair = Ndummy by definition,
              ! we already know N*(N-1) = 2*Ndummy
              if (i==nbins) print*, ij_type, Ndummy(ij_type)
              if (ij_type .ne. 5) then
                 ! the true pair correlation function is defined as
                 ! (A/N)*dN/(N-1)/(dA), where for us g is currently
                 ! dN, and dA = pi(r+dr)^2-pi(r)^2 = 2pi r dr+ pi dr*dr
                 g(i,ij_type) = g(i,ij_type) / (2*Ndummy(ij_type)/(box(1)*box(2)) ) / &
                      ( pi*(2*i+1)*dr*dr  ) 
              else
                 ! since for the border type, there is no double counting as
                 ! each side of the pair is distinct, N*(N-1) = Ndummy
                 g(i,ij_type) = g(i,ij_type) / ( Ndummy(ij_type) /&
                      (box(1)*box(2)) ) / ( pi*(2*i+1)*dr*dr  ) / 2
              endif
           endif
        enddo
        
     enddo
     

     do kx=1, nbinsq
        do ky=1, nbinsq
           sq(kx,ky,1) = qmin + (kx+0.5)*dq
           sq(kx,ky,2) = qmin + (ky+0.5)*dq
           ! not sure what's happening here, it looks s(q) is defined as
           ! s(q) = q_y + sum_{i=1}^{N-1}sum_{j=1}^{N-1} cos(q(r_i-r_j))
           ! (after some trig identities).
           sq(kx,ky,3) = sq(kx,ky,2) + sin_rq(kx,ky)*sin_rq(kx,ky) + cos_rq(kx,ky)*cos_rq(kx,ky)
        enddo
     enddo
     sq(:,:,3) = sq(:,:,3) / N 
     do i=2,5
        if (Ndummy(i)>0)  s_pr(i) = s_pr(i)/Ndummy(i) 
     enddo
     
     
  end if

end subroutine distance_matrix

subroutine get_cluster_MSD(ucoords0, ucoords1, box,  cl_i,&
     cl_s, N, D, Nc, MSD)
  !===========================================================
  ! Compute the averaged mean squared displacement of
  ! all clusters between two
  ! time points.
  !
  ! Returns the scalar mean squared displacement <R(t-t')^2>
  !
  ! cl_i is a list of all the atoms, each atom being labelled
  ! by the cluster it is in, with label 0 corresponding to
  ! the largest cluster and label Nc-1 corresponding to the
  ! smallest cluster.
  ! cl_s is the list of all cluster sizes, ordered from
  ! largest to smallest.
  ! N is the number of atoms in the simulation
  ! D is the dimension of the simulation
  ! Nc is the number of clusters
  !===========================================================
  implicit none
  integer, intent(in)   :: N, D, Nc
  real(8), intent(in)   :: ucoords0(N,D), ucoords1(N,D), box(D)
  integer, intent(in)   :: cl_i(N), cl_s(Nc)
  real(8), intent(out)  :: MSD 
  real(8)               :: com0(D), com1(D), delta, x(N), y(N)
  real(8)               :: tmp0(N,D), tmp1(N,D), r(D)
  integer               :: i, j, k
  
  MSD = 0
  ! loop over clusters
  do i=0, Nc-1
     ! count number of atoms in the cluster
     k = 0
     do j=1, N
        if (cl_i(j) == i) then
           k=k+1
           tmp0(k, :) = ucoords0(j,:)
           tmp1(k, :) = ucoords1(j,:)
        endif
        if (k == cl_s(i+1)) then
           exit
        endif
     enddo
     !Compute centre of masses of the clusters
     call get_com(tmp0, box, k, N, D, com0)
     call get_com(tmp1, box, k, N, D, com1)
     r = com1 - com0
     MSD = MSD + sum(r*r)
  enddo
  MSD = MSD / Nc
end subroutine get_cluster_MSD


subroutine get_cluster_omega(ucoords, v, box,  cl_i,  cl_s,&
     N, D, Nc, c_omega, c_v  )
  !===========================================================
  ! Compute angular momentum and linear momentum (speed) of
  ! each cluster.
  !
  ! Returns c_omega which is an array of length Nc storing
  ! angular momentum of each cluster, and c_v which is the
  ! same but for speed.
  !
  ! cl_i is a list of all the atoms, each atom being labelled
  ! by the cluster it is in, with label 0 corresponding to
  ! the largest cluster and label Nc-1 corresponding to the
  ! smallest cluster.
  ! cl_s is the list of all cluster sizes, ordered from
  ! largest to smallest.
  ! N is the number of atoms in the simulation
  ! D is the dimension of the simulation
  ! Nc is the number of clusters
  !===========================================================
  implicit none
  integer, intent(in)   :: N, D, Nc
  real(8), intent(in)   :: ucoords(N,D), box(D), v(N,D)
  integer, intent(in)   :: cl_i(N), cl_s(Nc)
  real(8), intent(out)  :: c_omega(Nc), c_v(Nc)
  real(8)               :: com(D), delta, x(N), y(N)
  real(8)               :: tmp(N,D), r(D), tmpv(N,D), com_v(D)
  
  integer            :: i, j, k, m
  
  c_omega = 0
  c_v = 0 
  ! loop over all clusters
  do i=0, Nc-1
     k = 0   ! count of atoms in cluster
     ! loop over all atoms
     do j=1, N
        if (cl_i(j) == i) then
           ! if jth atom is in cluster i
           k=k+1
           tmp(k, :) = ucoords(j,:)
           tmpv(k, :) = v(j,:)
        endif
        if (k == cl_s(i+1)) then
           ! if count of atoms in cluster i equals
           ! the number of atoms in cluster i,
           ! we're done.
           exit
        endif
     enddo

     c_omega(i+1) = 0
     com_v = 0

     ! compute centre of mass for cluster
     call get_com(tmp, box, k, N, D, com)
     ! loop over all atoms in the cluster
     do m=1, k
        ! compute distance of atom from centre of mass
        r = tmp(m,:) - com
        ! recast into periodic boundary condition coords
        call PBC(r, box)
        ! compute angular momentum contribution of this atom
        ! and add it to the cluster angular momentum.
        c_omega(i+1) = c_omega(i+1) + ( r(1) * tmpv(m,2) &
             - r(2) * tmpv(m,1) )
        ! compute velocity contribution of the particle and
        ! add it to the total velocity
        com_v = com_v + tmpv(m,:) 
     enddo
     !com_v = com_v / k
     ! compute speed of particle as centre of mass velocity squared
     c_v(i+1) = dsqrt(sum(com_v * com_v))
  enddo
end subroutine get_cluster_omega



subroutine get_overlap(ucoords0, ucoords1, psi6_re0,&
     psi6_im0,psi6_re1, psi6_im1, box, N, D, Q, g6t,&
     g6t_re, g6t_im)
  !===========================================================
  ! Compute the overlap autocorrelation function, defined as
  !
  !  Q(t,t') = 1/N sum_{i=1}^N exp(-[r_i(t+t')-r_i(t)]^2/a^2)
  !
  ! where r_i(t+t') corresponds to the coordinates
  ! ucoords1(i,:) and r_i(t) to ucoords0(i,:). Note that the
  ! coordinates are in REAL SPACE, NOT periodic boundary
  ! condition space.
  !
  ! Also compute the psi6 autocorrelation function,
  !
  !     1/N sum_{i=1}^N psi6(r_i(t+t'))psi6(r_i(t))*,
  !
  ! and return its real part g6t_re), imaginary part
  ! g6t_im, and modulus g6t.
  !
  ! N is number of atoms in simulation, D is dimension.
  ! ucoords0 is coordinates of one cluster, ucoords1 is other.
  ! box is size of periodic boundary conditions box.
  ! psi6 terms are real and imaginary parts of per atom hexagonal
  ! order parameter of clusters 1 and 2.
  !
  !===========================================================
  implicit none
  integer, intent(in)   :: N, D
  real(8), intent(in)   :: ucoords0(N,D), ucoords1(N,D), box(D)
  real(8), intent(in)   :: psi6_re0(N), psi6_im0(N)
  real(8), intent(in)   :: psi6_re1(N), psi6_im1(N)
  real(8), intent(out)  :: Q, g6t_re, g6t, g6t_im
  integer               :: i, j, k
  real(8)               :: r(D), rsq, a, asq, re, im
  
  a = 1.0
  asq = a * a
  Q = 0
  g6t_re = 0
  g6t_im = 0
  g6t = 0
  
  do i=1, N
     r = ucoords1(i,:) - ucoords0(i,:) 
     !call PBC(r, box)
     rsq = sum(r*r)
     Q = Q + exp(-rsq/asq)
     re = psi6_re0(i)*psi6_re1(i) + psi6_im0(i)*psi6_im1(i)
     im = - psi6_re0(i)*psi6_im1(i) + psi6_re1(i)*psi6_im0(i)
     g6t_re = g6t_re + re 
     g6t_im = g6t_im + im 
     g6t = g6t + dsqrt(re*re + im*im)
  enddo
  
  Q      = Q      / N
  g6t_re = g6t_re / N
  g6t_im = g6t_im / N
  g6t    = g6t    / N
  
end subroutine get_overlap


subroutine get_com(ucoords, box, Np, N, D, com)
  !===========================================================
  ! Compute centre of mass location for a given set of Np
  ! coordinates (from a subset of N total coordinates of
  ! dimension D).
  !===========================================================
  implicit none
  integer, intent(in)      :: N, D, Np
  real(8), intent(inout)   :: ucoords(N,D)
  real(8), intent(in)      :: box(D)
  real(8), intent(out)     :: com(D) 
  real(8)                  :: r(D)
  integer                  :: i
  do i=2, Np
     ! compute real space difference between two sets of coords
     r = ucoords(i, :) - ucoords(1, :)
     call PBC(r, box)
     ! change this to periodic boundary space difference
     ucoords(i, :) = ucoords(1, :) + r
  end do
  ! compute centre of mass coordinate of coords
  com = sum(ucoords(1:Np,:), DIM=1)/ Np
end subroutine get_com



subroutine dmatrix(coords, box, cutoff, MAXnb,&
     N, D,distance, distance_x,distance_y,num_nb,&
     list_nb)
  !===========================================================
  ! This function computes lots of things at the same time
  ! for efficiency.
  !
  ! Compute distance(N,N), which is the pairwise distance
  !  between all of the atoms.
  !
  ! Compute distance_x(N,N), distance_y(N,N) which are just
  !  x and y components of distance between pairs.
  !
  ! Compute list_nb(N,MAXnb) for which each atom tracks
  !  a list of its neighbours.
  !
  ! Compute a num_nb(N), a count for the number of neighbours
  !  each atom has.
  !
  ! N is number of atoms in the simulation.
  ! D is dimension.
  ! MAXnb is maximum number of neighbours expected (100 is
  !  a safe bet).
  ! coords(N,D) is the coordinates of all atoms.
  ! box(D) is the dimensions of the periodic box.
  ! cutoff is the distance at which two atoms are no
  !  longer neighbours.
  !===========================================================
  implicit none
  integer, intent(in)  :: N, D, MAXnb
  real(8), intent(in)  :: coords(N,D)
  real(8), intent(in)  :: box(D), cutoff
  real(8), intent(out) :: distance(N, N)
  real(8), intent(out) :: distance_x(N, N),distance_y(N, N)
  integer, intent(out) :: list_nb(N, MAXnb), num_nb(N)
  integer              :: i, j
  real(8)              :: rsq, dis, r(D), cutoff_sq
  
  distance = 0.0
  num_nb = 0
  cutoff_sq = cutoff *  cutoff 


  do i=1, N-1
     ! sum over atoms 1 to N-1     

     do j=i+1, N
        ! combined with above sum over 
        ! sum over each pair of atoms
        r = coords(i,:) - coords(j,:) 
        call PBC(r,box) 
        rsq = sum(r*r)
        dis = dsqrt(rsq)
        if (rsq <= cutoff_sq) then
           ! increase neighbour count for both i and j atoms
           num_nb(i) = num_nb(i) + 1
           num_nb(j) = num_nb(j) + 1
           ! track where neighbours are in list
           list_nb(i, num_nb(i)) = j
           list_nb(j, num_nb(j)) = i
        end if
        ! track neighbour matrix for distance
        distance(i,j) = dis
        distance(j,i) = dis
        ! track neighbour matrix for positions
        distance_x(i,j) = -r(1)
        distance_x(j,i) = r(1)      
        distance_y(i,j) = -r(2)
        distance_y(j,i) = r(2)        
        
     enddo
  enddo

end subroutine dmatrix


subroutine PBC(dd, box)
  !===========================================================
  ! Change coordinates from real space to periodic boundary
  ! conditions space.
  !
  ! dd is an array of length 2, which has the current (real
  ! space) coordinates of an atom.
  !
  ! box is the dimensions of the box (i.e. 100 by 100) in
  ! array form
  !===========================================================
  implicit none
  real(8), intent(in)    :: box(2)
  real(8), intent(inout) :: dd(2)
  integer i
  do i=1, 2
     do while (dd(i)>box(i)/2)
        dd(i) = dd(i) - box(i)
     end do
     do while (dd(i)<-box(i)/2)
        dd(i) = dd(i) + box(i)
     end do
  enddo
end subroutine PBC
