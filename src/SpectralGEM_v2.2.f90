
module write_module
	integer(4), parameter						:: lline=500
	character(lline)							:: wline
end module write_module

module input_module
	character(5)								:: ext
	integer(4)									:: min_cluster
	integer(4)									:: output
	integer(4)									:: l_id
	character(150)								:: data_dir
	integer(4)									:: l_dir
	character(50)								:: MMp_file
	character(50)								:: exclusion_file
end module input_module

module parameter_module
	integer(4), parameter						:: lid=35
	integer(4), parameter						:: min_sig_ev=3   !minimum of min_sig_ev-1 dimensions
	integer(4)									:: max_sig_ev  !maximum of max_sig_ev-1 dimensions
	integer(4)									:: n_snp_use	!number of snp to use in calculations
end module parameter_module	

module data_module
	use parameter_module
	integer(4)									:: n_ind	!total number of individuals
	integer(4)									:: n_snp	!number of SNP used
	type individual_data
		character(lid)							:: id		!indentification
		integer(4)								:: sex		!sex
		integer(4)								:: dx		!case control status
		integer(4)								:: grp		!additional grouping
		integer(4)								:: incl		!indicator whether individuals is included
		integer(4)								:: sub		!indicator whether individuals is included in the subset
		integer(4)								:: cluster	!indicator in which cluster an individual resides, -1 indicates that individual is not included
	end type individual_data
	type(individual_data), dimension(:), allocatable :: ind	!information on the individuals
	real(8), dimension(:,:), allocatable		:: MMp


end module data_module

module work_module
	use parameter_module
!	use dfport
!	real									:: tb, te							

	integer(4), dimension(:), allocatable		:: list
	character(lid), dimension(:), allocatable	:: id_list
	real(8), dimension(:), allocatable			:: eval
	real(8), dimension(:,:), allocatable		:: evec
	real(8), dimension(:,:), allocatable		:: dist
	real(8), dimension(:,:), allocatable		:: min_distance
	real(8)										:: n_eff
	integer(4)									:: n_cluster
	integer(4)									:: master_ev
	real(8), dimension(:), allocatable			:: sd_ev


end module work_module


subroutine genetic_distance_matching(input_file,matching_answer)
!this is the main program, it uses the idea of recalculating the eigenspace and distances after 
!each cluster that is formed.

!This is version 2.1 (5/18/2009)
!In this version there are a number of changes with the previous version
! 1) A bug was fixed that sometimes caused the program to crash
! 2) A new formula is used for determing the critical value of the eigengap.  This new
!    value now depends on the number of SNPs used for calculating the MM' matrix.  This allows
!    more accurate results in cases where there are few SNPs
! 3) The program will now use a minimum of 2 dimensions and a maximum of 50 dimensions.
!
! Bert Klei
! Computational Genetics
! Western Psychiatric Institute and Clinic
! University of Pittsburgh Medical Center
! kleil at upmc dot edu
!
!
 
use write_module
use data_module
use work_module
implicit none
!integer(4)										:: nsig_ev
!integer(4)										:: i_ind
character(100)                                          :: input_file
character(1)									:: matching_answer
character(1)									:: cluster_answer

n_eff=0
!call cpu_time(tb)

!read the input data
call read_input_data(input_file)

!read the MMp matrix
call read_MMp_matrix

write(wline,'(a)')																				;  call writeline
write(wline,'(a)')'**************************************************************************'  ;  call writeline
write(wline,'(a,i10)')'maximum number of significant dimensions :: ',max_sig_ev					;  call writeline
write(wline,'(a,i10)')'maximum SNP number used in calculations  :: ',n_snp_use					;  call writeline
write(wline,'(a)')'**************************************************************************'  ;  call writeline
write(wline,'(a)')																				;  call writeline

!transform MMp matrix
call transform_MMp_matrix

!check if there are individuals that need to be excluded
call exclude_individuals

!determine dimension of the problem
call determine_dimensions

write(wline,'(a)')	; call writeline
write(wline,'(a)')'============================================================================'	; call writeline
write(wline,'(a)')'For determinging ANCESTRY CLUSTERS and outliers enter C                     ' ; call writeline
write(wline,'(a)')'For determinging MATCHES for conditional logit association analysis enter M ' ; call writeline
write(wline,'(a)')'============================================================================'	; call writeline
!write(6,'(a)',advance='no')'enter C or M        :: '
!read(5,'(a)')matching_answer
write(wline,'(a)')	; call writeline
if(matching_answer == 'c' .or. matching_answer == 'C')then
	write(wline,'(a)')	; call writeline
	write(wline,'(a)')'CONTINUING WITH DETERMINING CLUSTERS and outliers'; call writeline
!	call cpu_time(tb)
	call determine_clusters


!	call cpu_time(te)
!	print *,'CPU time used for clustering :: ',te-tb

	call scale_eigenvectors

	call write_minimum_distance_data

!	call write_minimum_distance_R

!	write(wline,'(a)')	; call writeline
!	write(wline,'(a)')	; call writeline
!	write(wline,'(a)')'============================================================================'	; call writeline
!	write(wline,'(a)')'determine the cut-off to use through simulation '	; call writeline
!	write(wline,'(a,i5,a,i5,a,i4,a)')'use the command max.funNJW(',count(ind%dx == 2 .and. ind%incl == 1) &
!										,',',count(ind%dx == 1 .and. ind%incl == 1),',',master_ev,',  100)'; call writeline
!	write(wline,'(a)')'WARNING, the number of dimensions to use is 1 MORE than the SIGNIFICANT ones'; call writeline
!	write(wline,'(a)')'============================================================================'	; call writeline
!	write(wline,'(a)')	; call writeline
!	write(wline,'(a)')	; call writeline
	write(wline,'(a)')	; call writeline
else
	!create the files needed for the final matching

	call write_matching_files(master_ev)

!	call write_matching_R

	call significant_eigenvectors(master_ev,2)

	call write_eigenvalue_decomposition(master_ev)

endif
!call cpu_time(te)
!write(6,*)'cpu time need to complete the program :: ',te-tb
close(16)
call deallocate_arrays
return
end

subroutine deallocate_arrays
use data_module
use work_module

if(allocated(MMp)) deallocate(MMp)
if(allocated(ind)) deallocate(ind)
if(allocated(list)) deallocate(list)
if(allocated(id_list)) deallocate(id_list)
if(allocated(eval)) deallocate(eval)
if(allocated(evec)) deallocate(evec)
if(allocated(dist)) deallocate(dist)
if(allocated(min_distance)) deallocate(min_distance)
if(allocated(sd_ev)) deallocate(sd_ev)

return
end


subroutine scale_eigenvectors
!this subroutine scales the eigenvectors so that the final outliers can be identified
use write_module
use data_module
use work_module
use input_module
implicit none
integer(4)					:: m_ind, i_ind,l_ind
integer(4)					:: n_ev, i_ev
integer(4)					:: i_cluster
integer(4)					:: n_star
integer(4), dimension(:), allocatable		:: cluster_size
real(8), dimension(:,:), allocatable			:: mean_ev

write(wline,'(a)')	; call writeline
write(wline,'(a)')'RESCALE THE EIGENVECTORS BASED ON THE CLUSTERS'	; call writeline

allocate(cluster_size(1:n_cluster)); cluster_size=0
do i_ind=1,n_ind
	if(ind(i_ind)%cluster > 0)then
		cluster_size(ind(i_ind)%cluster)=cluster_size(ind(i_ind)%cluster)+1
	endif
enddo
write(wline,'(a,i6)')'total number of clusters                      :: ',count(cluster_size > 0)	; call writeline
write(wline,'(a,i6)')'number of HOMOGENEOUS clusters                :: ',count(cluster_size >= min_cluster)	; call writeline
write(wline,'(a,i6)')'number of       SMALL clusters                :: ',count(cluster_size <  min_cluster)	; call writeline
write(wline,'(a)')	; call writeline
write(wline,'(a,i6)')'total number of individuals                   :: ',count(ind%incl > 0)	; call writeline
write(wline,'(a,i6)')'number of individuals in HOMOGENEOUS clusters :: ',sum(cluster_size, MASK= cluster_size >= min_cluster)	
						 call writeline
write(wline,'(a,i6)')'number of individuals in SMALL clusters       :: ',sum(cluster_size, MASK= cluster_size < min_cluster)	
						call writeline


!make a list of individuals to be included
ind%sub=0; m_ind=0
do i_ind=1,n_ind
	if(ind(i_ind)%incl == 1)then
		ind(i_ind)%sub=1
		m_ind=m_ind+1
	endif
enddo
write(wline,'(a,i6)')'number of individuals for which eigenvectors are to be rescaled :: ',m_ind	; call writeline
!determine the significant eigenvectors
n_ev=master_ev
call significant_eigenvectors(n_ev,2)

write(wline,'(a)')	; call writeline
write(wline,'(a)')'write the file with cluster assignments'	;call writeline

open(100,file=data_dir(1:l_dir)//'clusters_'//ext//'.txt')
l_ind=0

do i_ind=1,n_ind
	if(ind(i_ind)%cluster > -1)then
		l_ind=l_ind+1
		write(100,'(a,3i3,i6,100f10.6)')ind(i_ind)%id(1:l_id),ind(i_ind)%sex,ind(i_ind)%dx,ind(i_ind)%grp,ind(i_ind)%cluster, &
															evec(l_ind,1:n_ev)
	endif
enddo
close(100)
write(wline,'(a)')'cluster assignments can be found in file (format: id, sex, dx, cluster, eigenvectors) :: '	; call writeline
write(wline,'(8a)')'        ',data_dir(1:l_dir),'clusters_',ext,'.txt'	; call writeline

call write_eigenvalue_decomposition(n_ev)

write(wline,'(a)')	; call writeline
write(wline,'(a)')'determine the average eigenvector loading for each cluster'	; call writeline

!determine the average eigenvector load
allocate(mean_ev(n_cluster,n_ev)); mean_ev=0
!find the sum for each cluster
do i_ind=1,m_ind
	do i_ev=1,n_ev
		mean_ev(ind(list(i_ind))%cluster,i_ev)=mean_ev(ind(list(i_ind))%cluster,i_ev)+evec(i_ind,i_ev)
	enddo
enddo

!now calculate the average
mean_ev(:,1)=0
do i_cluster=1,n_cluster
	mean_ev(i_cluster,:)=mean_ev(i_cluster,:)/cluster_size(i_cluster)
enddo


!determine the standard deviation
write(wline,'(a)')'determine the standard deviation for the eigenvecor loading in the large clusters'	; call writeline
allocate(sd_ev(n_ev)); sd_ev=0
do i_ind=1,m_ind
	if(cluster_size(ind(list(i_ind))%cluster) >= min_cluster)then
		do i_ev=1,n_ev
			sd_ev(i_ev)=sd_ev(i_ev)+(evec(i_ind,i_ev)-mean_ev(ind(list(i_ind))%cluster,i_ev))**2
		enddo
	endif
enddo

!adjust for the individuals that are not used
n_star=sum(cluster_size, MASK= cluster_size >= min_cluster)-count(cluster_size >= min_cluster)
write(wline,'(a)')
write(wline,'(a)')'adjust the number of individuals for small groups and number of homogenous groups'	; call writeline
write(wline,'(a,i10)')'number of individuals                    :: ',m_ind	; call writeline
write(wline,'(a,i10)')'number of indivdiuals in large clusters  :: ',sum(cluster_size, MASK= cluster_size >= min_cluster)	
																												call writeline
write(wline,'(a,i10)')'number of large clusters                 :: ',count(cluster_size >= min_cluster)	; call writeline
write(wline,'(a,i10)')'adjusted number of ind in large clusters :: ',n_star	; call writeline


sd_ev=sqrt(m_ind*sd_ev/n_star)
sd_ev(1)=1.d0
write(wline,'(a)')	; call writeline

!determine distance matrix
call calculate_distance_NJW(n_ev,2,2)

return
end

subroutine determine_clusters
!this subroutine assigns the clusters to the individuals
!algorithm finds the largest cluster without any significant eigenvalues based on the distance matrix
!after a group of individuals have been assigned the correlation matrix is decomposed again with these
!individuals removed.
use write_module
use data_module
use work_module
use input_module

implicit none
integer(4)								:: n_ev
integer(4)								:: i_ind, i
integer(4)								:: i_cluster, j_cluster
integer(4), dimension(:), allocatable	:: master_list
integer(4)								:: m_ind
integer(4), dimension(:), allocatable	:: old_cluster
integer(4), dimension(1)				:: c_clus

!clustering subroutine variables
integer(4)								:: imeth
integer(4)								:: idist
!cluster assignment
integer(4), dimension(:), allocatable	:: nclus
integer(4), dimension(:), allocatable	:: iclus
integer(4), dimension(:), allocatable	:: ix
integer(4), dimension(:), allocatable	:: iy

integer(4)								:: len
integer(4)								:: iopt
integer(4), dimension(:), allocatable	:: ia
integer(4), dimension(:), allocatable	:: ib
real(8), dimension(:), allocatable		:: crit
real(8), dimension(:), allocatable		:: membr
integer(4), dimension(:), allocatable	:: nn
real(8), dimension(:), allocatable		:: disnn
logical, dimension(:), allocatable		:: flag
real(8), dimension(:), allocatable		:: diss

integer(4), dimension(:,:), allocatable	:: iclass
integer(4), dimension(:), allocatable	:: hvals
integer(4), dimension(:), allocatable	:: iorder
real(8), dimension(:), allocatable		:: critval
integer(4), dimension(:), allocatable	:: height
integer(4)								:: i_level
real(8), dimension(:,:), allocatable	:: norm_evec
real(8)									:: norm


n_cluster=0
write(wline,'(a)')	; call writeline
write(wline,'(a)')'ASSIGN THE CLUSTERS TO THE INDIVIDUALS'	; call writeline
write(wline,'(a,i6)')   'minimum cluster size required                 :: ',min_cluster	; call writeline
write(wline,'(a,i6)')   'number of individuals to be clustered         :: ',count(ind%incl == 1)	; call writeline


m_ind=size(list)
write(wline,'(a,i6)')   'base number of significant eigenvalues        :: ',master_ev-1	; call writeline
n_ev=master_ev
do while(count(ind%cluster == 0) >= min_cluster)
	if(count(ind%cluster == 0) /= count(ind%incl == 1))then
		!flag the individuals to be included
		ind%sub=0; m_ind=0
		do i_ind=1,n_ind
			if(ind(i_ind)%incl == 1 .and. ind(i_ind)%cluster == 0)then
				ind(i_ind)%sub=1
				m_ind=m_ind+1
			endif
		enddo
		if(output == 2)then
			write(wline,'(a,i6)')'remaining number of individuals to be clustered :: ',m_ind	; call writeline
		endif
			
		!determine eigenvalue decomposition
		call determine_eigenvalues(output)
		!determine the significant eigenvalues
		call significant_eigenvalues(eval,m_ind,n_ev,output)

		if(n_ev == 1)exit
		!determine significant eigenvectors
		call significant_eigenvectors(n_ev,output)


	endif
	if(n_ev == 1)then
		write(wline,'(a)')'remaining group of individuals is HOMOGENEOUS'
		exit
	endif
	if(count(ind%sub == 1) < min_cluster)then
		write(wline,'(a)')'remaining group of individuals is TOO SMALL'
		exit
	endif

	!store the master list of individuals
	if(allocated(master_list))deallocate(master_list)
	allocate(master_list(size(list)));master_list=list


	!find the cluster tree for these individuals
	imeth=4
	idist=0
	m_ind=count(ind%sub == 1)
	len=(m_ind*(m_ind-1))/2

	if(allocated(ia))deallocate(ia)
	allocate(ia(m_ind))
	if(allocated(ib))deallocate(ib)
	allocate(ib(m_ind))
	if(allocated(crit))deallocate(crit)
	allocate(crit(m_ind))
	if(allocated(membr))deallocate(membr)
	allocate(membr(m_ind))
	if(allocated(nn))deallocate(nn)
	allocate(nn(m_ind))
	if(allocated(disnn))deallocate(disnn)
	allocate(disnn(m_ind))
	if(allocated(flag))deallocate(flag)
	allocate(flag(m_ind))
	if(allocated(diss))deallocate(diss)
	allocate(diss(len))
	iopt=1

	!norm the eigenvectors so it can be used for the NJW distance
	allocate(norm_evec(m_ind,size(eval)))
	do i_ind=1,m_ind
		norm=sqrt(sum(evec(i_ind,:)**2))
		norm_evec(i_ind,:)=evec(i_ind,:)/norm
	enddo
	call hc(m_ind,n_ev,len,iopt,norm_evec(:,1:size(eval)),ia,ib,crit, &
			membr,nn,disnn,flag,diss)

	deallocate(norm_evec)

	!now check the actual clusters
	if(allocated(old_cluster))deallocate(old_cluster)
	allocate(old_cluster(size(master_list)));old_cluster=0
	do i_cluster=2,size(master_list)
		!determine the cluster membership based on i_cluster clusters
		if(allocated(nclus))deallocate(nclus)
		if(allocated(iclus))deallocate(iclus)

		allocate(nclus(i_cluster),iclus(size(master_list)))

		!assign the individuals to i_cluster clusters
		!allocate space needed to hcass
!			call cpu_time(tb)
		i_level=i_cluster
		if(allocated(iclass))deallocate(iclass)
		allocate(iclass(m_ind,i_level))
		if(allocated(hvals))deallocate(hvals)
		allocate(hvals(i_level))
		if(allocated(iorder))deallocate(iorder)
		allocate(iorder(i_level))
		if(allocated(critval))deallocate(critval)
		allocate(critval(i_level))
		if(allocated(height))deallocate(height)
		allocate(height(i_level))
	
		iclass=0
		call hcass(m_ind,ia,ib,crit,i_level,iclass,hvals,iorder,critval,height)
	
		iclus=iclass(:,i_cluster-1)
		do i=1,i_cluster
			nclus(i)=count(iclass(:,i_cluster-1) == i)
		enddo



		!sort the clusters from largest to smallest
		if(allocated(ix))deallocate(ix)
		if(allocated(iy))deallocate(iy)
		allocate(ix(i_cluster),iy(i_cluster))
		ix=nclus
		do i=1,i_cluster;iy(i)=i;enddo
		call isort(ix,iy,i_cluster,-2)
		
		!determine the subset of individuals in this cluster
		do j_cluster=1,i_cluster
			!check if the cluster composition has changed
			c_clus=old_cluster(minloc(iclus, MASK= iclus == iy(j_cluster)))
			if(count(iclus == iy(j_cluster)) == count(old_cluster == c_clus(1)))then
				if(output == 2)then
					write(wline,'(a)')'cluster not divided'	; call writeline
					write(wline,'(a,i6,a,i6)')'previous count :: ',count(iclus == iy(j_cluster)),' new count :: ', &
													count(old_cluster == c_clus(1)) ; call writeline
				endif
				cycle
			endif

			ind%sub=0
			do i_ind=1,size(master_list)
				if(iclus(i_ind) == iy(j_cluster))then
					ind(master_list(i_ind))%sub=1
				endif 
			enddo
			if(count(ind%sub == 1) < min_cluster)then
				if(output == 2)then
					write(wline,'(a)')	; call writeline
				endif
				write(wline,'(a,i6,a,i6,a,i6)')'      Small cluster :: ',n_cluster+1,' # individuals :: ', &
					count(ind%sub == 1),' remaining :: ',count(ind%cluster == 0 .and. ind%sub == 0)	
					call writeline
				exit
			endif

			!determine eigenvalue decomposition for this subset
			call determine_eigenvalues(output)

			!determine the significant eigenvalues
			call significant_eigenvalues(eval,size(eval),n_ev,output)

			if(n_ev == 1)then
				if(output == 2)then
					write(wline,'(a)')	; call writeline
				endif
				write(wline,'(a,i6,a,i6,a,i6)')'Homogeneous cluster :: ',n_cluster+1,' # individuals :: ', &
					count(ind%sub == 1),' remaining :: ',count(ind%cluster == 0 .and. ind%sub == 0)
					call writeline
				exit
			endif
		enddo

		!assign cluster number to the individuals
		if(n_ev == 1 .or. count(ind%sub == 1) < min_cluster)then
			n_cluster=n_cluster+1
			do i_ind=1,n_ind
				if(ind(i_ind)%sub == 1)ind(i_ind)%cluster=n_cluster
			enddo
			if(output == 2)then
				write(wline,'(a,i6)')'number of individuals to be clustered  :: ', &
					count(ind%cluster == 0)	; call writeline
			endif
			exit
		endif
		old_cluster=iclus
	enddo

enddo


!assign clusters to any remaining individuals
if(count(ind%cluster == 0) > 0)then
	if(output == 2)then
		write(wline,'(a)')	; call writeline
	endif
	if(n_ev == 1 .and. count(ind%cluster == 0) >= min_cluster )then
		write(wline,'(a,i6,a,i6,a,i6)')'Homogeneous cluster :: ',n_cluster+1,' # individuals :: ', &
			count(ind%sub == 1),' remaining :: ',count(ind%cluster == 0 .and. ind%sub == 0)
			call writeline
	else
		write(wline,'(a,i6,a,i6,a,i6)')'      Small cluster :: ',n_cluster+1,' # individuals :: ', &
			count(ind%cluster == 0 .and. ind%sub == 0),' remaining :: ',0	
			call writeline
	endif
	n_cluster=n_cluster+1
	do i_ind=1,n_ind
		if(ind(i_ind)%cluster == 0)then
			ind(i_ind)%cluster=n_cluster
		endif
	enddo
endif

write(wline,'(a)')	; call writeline
write(wline,'(a)')'CLUSTERING COMPLETED'	; call writeline

if(allocated(nclus))deallocate(nclus)
if(allocated(iclus))deallocate(iclus)
if(allocated(ix))deallocate(ix)
if(allocated(iy))deallocate(iy)

return
end

subroutine calculate_distance_NJW(n_sig,opt_out,match_step)
!This subroutine calculates the distances using the NJW algorithm
use data_module
use input_module
use write_module
use work_module
implicit none
integer(4)											:: i_ind, j_ind
integer(4)											:: n_sig, i_sig
integer(4)											:: m_ind
integer(4)											:: opt_out
integer(4)											:: n_ev
integer(4)											:: match_step

real(8)												:: norm
real(8),dimension(:,:),allocatable					:: norm_evec


if(opt_out == 2)then
	write(wline,'(a)')	; call writeline
	write(wline,'(a)')'CALCULATE THE DISTANCE MATRIX BASED ON THE SIGNIFICANT EIGENVECTORS USING THE NJW ALGORITHM'	; call writeline
endif

!determine the number of individuals
m_ind=count(ind%sub == 1)

!determine the size of the problem
n_ev=size(eval)

!now norm the eigenvectors
allocate(norm_evec(m_ind,n_ev))
do i_ind=1,m_ind
	norm=sqrt(sum(evec(i_ind,:)**2))
	norm_evec(i_ind,:)=evec(i_ind,:)/norm
enddo

if(.not. allocated(sd_ev))then
	allocate(sd_ev(n_ev)); sd_ev=1
endif


!scale the eigenvectors
do i_ind=1,m_ind
	do i_sig=1,n_sig
		norm_evec(i_ind,i_sig)=norm_evec(i_ind,i_sig)/sd_ev(i_sig)
	enddo
enddo


if(allocated(dist))deallocate(dist)
allocate(dist(m_ind,m_ind)); dist=0.d0

!now calculate the distances
do i_ind=1,m_ind
	do j_ind=1,i_ind
		do i_sig=1,n_sig
			dist(i_ind,j_ind)=dist(i_ind,j_ind)+(norm_evec(i_ind,i_sig)-norm_evec(j_ind,i_sig))**2
		enddo
		dist(i_ind,j_ind)=sqrt(dist(i_ind,j_ind))
		dist(j_ind,i_ind)=dist(i_ind,j_ind)
	enddo
enddo


!determine the minimum distances between cases and controls
if(allocated(min_distance))deallocate(min_distance)
allocate(min_distance(n_ind,2)); min_distance=99999
do i_ind=1,m_ind
	do j_ind=1,m_ind
		if(ind(list(i_ind))%dx /= ind(list(j_ind))%dx)then
			min_distance(list(i_ind),1)=min(min_distance(list(i_ind),1),dist(i_ind,j_ind))
		elseif(i_ind /= j_ind)then
			min_distance(list(i_ind),2)=min(min_distance(list(i_ind),2),dist(i_ind,j_ind))
		endif
	enddo
enddo

if(opt_out == 2)then
	write(wline,'(a)')	; call writeline
	write(wline,'(a,f10.4)')'maximum distance between a case    and the closest control :: ', &
							maxval(min_distance(:,1), MASK = ind%dx == 2 .and. ind%sub == 1)	; call writeline
	write(wline,'(a,f10.4)')'maximum distance between a control and the closest case    :: ', &
							maxval(min_distance(:,1), MASK = ind%dx == 1 .and. ind%sub == 1)	; call writeline
	write(wline,'(a,f10.4)')'maximum distance between a case    and the closest case    :: ', &
							maxval(min_distance(:,2), MASK = ind%dx == 2 .and. ind%sub == 1)	; call writeline
	write(wline,'(a,f10.4)')'maximum distance between a control and the closest control :: ', &
							maxval(min_distance(:,2), MASK = ind%dx == 1 .and. ind%sub == 1)	; call writeline
endif

if(match_step == 1)then
	evec=norm_evec
endif

!deallocate space used
deallocate(norm_evec)


return
end

subroutine significant_eigenvalues(ev,m,n_sig,opt_out)
!this subroutine determines which of the eigenvalues are significant based on the eigengap
!heuristic -0.00016 + 2.7/m + 2.3/n_snp
!eigenvalues are in eval
!m is the number of eigenvalues
!n_sig is the number of significant eigenvalues 
use data_module
use input_module
use write_module
use work_module
implicit none
integer(4)						:: m
real(8), dimension(1:m)			:: ev
integer(4)						:: n_sig
integer(4)						:: opt_out
real(8)							:: crit

integer(4)						:: i

real(8), dimension(:), allocatable :: eigen_gap

!find the signficance level
crit=-0.00016 + 2.7/m + 2.3/n_snp

if(opt_out == 2)then
	write(wline,'(a)')	; call writeline
	write(wline,'(a)')'DETERMINE THE SIGNIFICANT NUMBER OF EIGENVALUES'	; call writeline
	write(wline,'(a)')'USING THE EIGENGAP HEURISTIC'	; call writeline
	write(wline,'(a,f10.6)')'critical value -0.00016 + 2.7/m + 2.3/n_snp :: ',crit	; call writeline
endif

allocate(eigen_gap(m-1))

!eigenvalues are listed from smallest to biggest
!skip over eigenvalues greater than 1
do i=1,m
	if(ev(i) > 1.d0)then
		ev(i)=1.d0
	endif
enddo

!determine the gaps
if(opt_out == 2)then
	open(120,file=data_dir(1:l_dir)//'eigen_gap_'//ext//'.txt')
endif
do i=1,m-1
	eigen_gap(i)=ev(i+1)-ev(i)
	if(opt_out == 2)then
		write(120,'(i5,2es20.10)')i,eigen_gap(i)
	endif
enddo
if(opt_out == 2)then
	close(120)
endif
!print '(6f10.6)',eigen_gap
!pause

n_sig=min(max_sig_ev,m-1)
do i=1,min(max_sig_ev,m-1)
	if(eigen_gap(i) > crit)then
		n_sig=i
	endif
enddo

deallocate(eigen_gap)
if(opt_out == 2)then
	write(wline,'(a)')	; call writeline
	write(wline,'(a,i6)')'NUMBER OF SIGNIFICANT EIGENVALUES :: ',n_sig-1	; call writeline
endif
!pause

return
end



subroutine significant_eigenvectors(n_sig,opt_out)
!this subroutine determines the eigenvectors for the n_sig significant eigenvalues
use data_module
use write_module
use input_module
use work_module
implicit none
integer(4)								:: n_sig
integer(4)								:: opt_out
integer(4)								:: m_ind, i_ind, j_ind

real(8), dimension(:,:), allocatable	:: MMp_sub
real(8), dimension(:), allocatable		:: DEG
integer(4)								:: m, info, lwork, liwork
integer(4), dimension(:), allocatable	:: isuppz
real(8), dimension(:), allocatable		:: work
integer(4), dimension(:), allocatable	:: iwork


if(opt_out == 2)then
	write(wline,'(a)')	; call writeline
	write(wline,'(a)')   'DETERMINE THE SIGNIFICANT EIGENVECTORS'	; call writeline
endif

m_ind=count(ind%sub == 1)
if(m_ind /= size(list))then
	if(allocated(list))deallocate(list)
	allocate(list(m_ind))
	if(allocated(id_list))deallocate(id_list)
	allocate(id_list(m_ind))
	m_ind=0
	do i_ind=1,n_ind
		if(ind(i_ind)%sub == 1)then
			m_ind=m_ind+1
			list(m_ind)=i_ind
			id_list(m_ind)=ind(i_ind)%id
		endif
	enddo
	if(opt_out == 2)then
		write(wline,'(a,i6)')'list of individuals created    :: ',m_ind	; call writeline
	endif
endif
if(opt_out == 2)then
	write(wline,'(a,i6)')'number of individuals included     :: ',m_ind ; call writeline
endif

!create the submatrix from the correlation matrix
allocate(MMp_sub(m_ind,m_ind))
if(allocated(evec))deallocate(evec)
allocate(evec(m_ind,n_sig))
if(allocated(eval))deallocate(eval)
allocate(eval(n_sig))

!create the submatrix from the correlation matrix
do i_ind=1,m_ind
	do j_ind=1,m_ind
		MMp_sub(i_ind,j_ind)=MMp(list(i_ind),list(j_ind))
	enddo
enddo
if(opt_out == 2)then
	write(wline,'(a)')'sub-matrix of MMp created '	; call writeline
endif

!create the degree vector
if(opt_out == 2)then
	write(wline,'(a)')	; call writeline
	write(wline,'(a)')'create the deg=1/sqrt(DEG) vector '	; call writeline
endif	

allocate(DEG(m_ind))
do i_ind=1,m_ind
	DEG(i_ind)=sqrt(1.d0/sum(MMp_sub(i_ind,:)))
enddo

!create lambda
if(opt_out == 2)then
	write(wline,'(a)')	; call writeline
	write(wline,'(a)')'create NORMALIZED LAPLACIAN=(I-deg_i*W(ij)*deg_j) '	; call writeline
endif	

do i_ind=1,m_ind
	do j_ind=1,i_ind
		MMp_sub(i_ind,j_ind)=-DEG(i_ind)*MMp_sub(i_ind,j_ind)*DEG(j_ind)
		MMp_sub(j_ind,i_ind)=MMp_sub(i_ind,j_ind)
	enddo
	MMp_sub(i_ind,i_ind)=MMp_sub(i_ind,i_ind)+1
enddo

if(opt_out == 2)then
	write(wline,'(a)')	; call writeline
	write(wline,'(a,i3,a)')'determine the ',n_sig-1,' significant eigenvectors of Normalized Laplacian'	; call writeline
endif	

!call cpu_time(tb) 
allocate(isuppz(2*n_sig)) !; isuppz=0
lwork=26*m_ind
allocate(work(lwork)) !; work=0
liwork=10*m_ind
allocate(iwork(liwork)) !; iwork=0

call dsyevr('V','I','L',m_ind,MMp_sub,m_ind,0,0,1,n_sig,-1.d0,m, &
			eval,evec,m_ind,isuppz,work,lwork,iwork,liwork,info)
!call cpu_time(te)
!print *,'lapack ',te-tb

!deallccate space
deallocate(isuppz,work,iwork)
deallocate(MMp_sub)
deallocate(DEG)


return
end

subroutine determine_eigenvalues(opt_out)
!this subroutine determines the dimension of the original problem
use data_module
use write_module
use input_module
use work_module
implicit none
integer(4)								:: opt_out
integer(4)								:: i_ind, j_ind
integer(4)								:: m_ind
real(8), dimension(:,:), allocatable	:: MMp_sub
real(8), dimension(:), allocatable		:: DEG
real(8), dimension(:), allocatable		:: work
integer(4)								:: lwork,info

if(opt_out == 2)then
	write(wline,'(a)')	; call writeline
	write(wline,'(a)')'DETERMINE THE EIGENVALUES'	; call writeline
endif

m_ind=count(ind%sub == 1)
if(opt_out == 2)then
	write(wline,'(a,i6)')'number of individuals included    :: ',m_ind	; call writeline
endif

if(allocated(list))deallocate(list)
allocate(list(m_ind))
if(allocated(id_list))deallocate(id_list)
allocate(id_list(m_ind))
m_ind=0
do i_ind=1,n_ind
	if(ind(i_ind)%sub == 1)then
		m_ind=m_ind+1
		list(m_ind)=i_ind
		id_list(m_ind)=ind(i_ind)%id
	endif
enddo
if(opt_out == 2)then
	write(wline,'(a,i6)')'list of individuals created    :: ',m_ind	; call writeline
endif

!create the submatrix from the correlation matrix
allocate(MMp_sub(m_ind,m_ind))
do i_ind=1,m_ind
	do j_ind=1,m_ind
		MMp_sub(i_ind,j_ind)=MMp(list(i_ind),list(j_ind))
	enddo
enddo
if(opt_out == 2)then
	write(wline,'(a)')'sub-matrix of MMp created '	; call writeline
endif

!create the degree vector
if(opt_out == 2)then
	write(wline,'(a)')	; call writeline
	write(wline,'(a)')'create the deg=1/sqrt(DEG) vector '	; call writeline
endif	

allocate(DEG(m_ind))
do i_ind=1,m_ind
	DEG(i_ind)=sqrt(1.d0/sum(MMp_sub(i_ind,:)))
enddo

!create lambda
if(opt_out == 2)then
	write(wline,'(a)')	; call writeline
	write(wline,'(a)')'create NORMALIZED LAPLACIAN = (I-deg_i*W(ij)*deg_j) '	; call writeline
endif	

do i_ind=1,m_ind
	do j_ind=1,i_ind
		MMp_sub(i_ind,j_ind)=-DEG(i_ind)*MMp_sub(i_ind,j_ind)*DEG(j_ind)
		MMp_sub(j_ind,i_ind)=MMp_sub(i_ind,j_ind)
	enddo
	MMp_sub(i_ind,i_ind)=MMp_sub(i_ind,i_ind)+1
enddo
	

!set up space needed
if(allocated(eval))deallocate(eval)
allocate(eval(m_ind)); eval=0

if(opt_out == 2)then
	write(wline,'(a)')	; call writeline
	write(wline,'(a)')'determine eigenvalues of Normalized Laplacian '	; call writeline
endif	

!call cpu_time(tb)
lwork=3*m_ind-1
allocate(work(lwork))

call dsyev('N','L',m_ind,MMp_sub,m_ind,eval,work,lwork,info)

!call cpu_time(te)
!print *,'lapack time :: ',te-tb

!deallocate the space used
deallocate(work)
deallocate(MMp_sub)
deallocate(DEG)

return
end

subroutine determine_dimensions
!this subroutine important dimensions of the problem
use write_module
use data_module
use input_module
use work_module
implicit none
integer(4)				:: n_ev
integer(4)				:: m_ind


write(wline,'(a)')	; call writeline
write(wline,'(a)')'DETERMINE THE DIMENSION OF THE PROBLEM'	; call writeline
ind%sub=ind%incl

!determine the set of eigenvalues for the complete dataset
call determine_eigenvalues(2)

!determine the number of significant eigenvalues
m_ind=count(ind%sub == 1)
call significant_eigenvalues(eval(1:m_ind),m_ind,n_ev,2)

master_ev=n_ev
if(master_ev /= min(max(min_sig_ev,n_ev),max_sig_ev))then
	master_ev=min(max(min_sig_ev,n_ev),max_sig_ev)
	write(wline,'(a)')												; call writeline
	write(wline,'(a)')'********************************************'; call writeline
	write(wline,'(a,i6)')'NUMBER OF EIGENVALUES USED FOR CALCULATIONS :: ',master_ev-1; call writeline
	write(wline,'(a)')'********************************************'; call writeline
endif

!determine the set of eigenvectors
call significant_eigenvectors(master_ev,2)

return
end


subroutine exclude_individuals
!this subroutine checks for the existence of a list of individuals that need to be 
!excluded from the analysis.  It marks these individuals with a 0 in the inclusion field of the
!individual information.
use data_module
use write_module
use input_module
implicit none
logical									:: file_exist
integer(4)								:: eof
integer(4)								:: n_rec
character(lid)							:: id
integer(4)								:: n_grp, i_grp
integer(4)								:: i_ind

write(wline,'(a)')	; call writeline
write(wline,'(a)')'CHECK IF ANY INDIVIDUALS NEED TO BE EXCLUDED FROM THE ANALYSIS'	; call writeline
write(wline,'(a)')'these need to be listed in the file :: '	; call writeline
write(wline,'(5x,a)')exclusion_file		; call writeline

ind%incl=1; ind%cluster=0
inquire(file=exclusion_file,exist=file_exist)
if(file_exist)then
	open(1,file=exclusion_file,status='old',action='read')
	eof=0; n_rec=0
	do
		read(1,*,iostat=eof)id
		if(eof < 0)exit
		n_rec=n_rec+1
		do i_ind=1,n_ind
			if(ind(i_ind)%id(1:l_id) == id(1:l_id))exit
		enddo
		if(i_ind <= n_ind)then
			ind(i_ind)%incl=0
			ind(i_ind)%cluster=-1
		else
			write(wline,'(a,a)')'Individual in exclusion list not found in MMp :: ',id(1:l_id); call writeline
		endif
	enddo
	write(wline,'(a,i6)')'Number of individuals excluded :: ',n_rec; call writeline

	close(1)
else
	write(wline,'(a)')	; call writeline
	write(wline,'(a)')'No individuals were excluded '	; call writeline
endif
write(wline,'(a)')	; call writeline
write(wline,'(a,i6)')'Number of individuals in the analysis      :: ',count(ind%incl == 1)	; call writeline 
write(wline,'(a,i6)')'Number of individuals excluded             :: ',count(ind%cluster < 0)	; call writeline 

write(wline,'(a)')	; call writeline
write(wline,'(a)')'COUNTS AFTER EXCLUSION OF INDIVIDUALS'	; call writeline
write(wline,'(a,i6)')'number of    males     :: ',count(ind%sex == 1 .and. ind%incl == 1)	; call writeline
write(wline,'(a,i6)')'number of  females     :: ',count(ind%sex == 2 .and. ind%incl == 1)	; call writeline
write(wline,'(a,i6)')'number of    cases     :: ',count(ind%dx == 2 .and. ind%incl == 1)	; call writeline
write(wline,'(a,i6)')'number of controls     :: ',count(ind%dx == 1 .and. ind%incl == 1)	; call writeline
n_grp =maxval(ind%grp)
write(wline,'(a)')	; call writeline
do i_grp=1,n_grp
	write(wline,'(a,i3,a,i6)')'number in group :: ',i_grp,' :: ',count(ind%grp == i_grp .and. ind%incl == 1)	; call writeline
enddo
write(wline,'(a)')	; call writeline
write(wline,'(a,120i6)')      '                      group ::        ',(i_grp,i_grp=1,n_grp)	; call writeline
write(wline,'(a,i6,1x,120i6)')'number of    male cases     :: ',count(ind%sex == 1 .and. ind%dx == 2 .and. ind%incl == 1), &
			(count(ind%sex == 1 .and. ind%dx == 2 .and. ind%grp == i_grp .and. ind%incl == 1),i_grp=1,n_grp)	; call writeline
write(wline,'(a,i6,1x,120i6)')'number of  female cases     :: ',count(ind%sex == 2 .and. ind%dx == 2 .and. ind%incl == 1), &
			(count(ind%sex == 2 .and. ind%dx == 2 .and. ind%grp == i_grp .and. ind%incl == 1),i_grp=1,n_grp)	; call writeline
write(wline,'(a,i6,1x,120i6)')'number of    male controls  :: ',count(ind%sex == 1 .and. ind%dx == 1 .and. ind%incl == 1), &
			(count(ind%sex == 1 .and. ind%dx == 1 .and. ind%grp == i_grp .and. ind%incl == 1),i_grp=1,n_grp)	; call writeline
write(wline,'(a,i6,1x,120i6)')'number of  female controls  :: ',count(ind%sex == 2 .and. ind%dx == 1 .and. ind%incl == 1), &
			(count(ind%sex == 2 .and. ind%dx == 1 .and. ind%grp == i_grp .and. ind%incl == 1),i_grp=1,n_grp)	; call writeline


return
end

subroutine transform_MMp_matrix
use data_module
use write_module
use input_module
implicit none
integer(4)							:: i_ind
integer(4)							:: j_ind

write(wline,'(a)')		; call writeline
write(wline,'(a)')'CONSTRUCT WEIGHT MATRIX FROM MMp'  ; call writeline

do i_ind=1,n_ind
	do j_ind=1,i_ind
		if(MMp(i_ind,j_ind) > 1.d-16)then
			MMp(i_ind,j_ind)=sqrt(MMp(i_ind,j_ind))
		else
			MMp(i_ind,j_ind)=0
		endif
		MMp(j_ind,i_ind)=MMp(i_ind,j_ind)
	enddo
enddo
write(wline,'(a)')'TRANSFORMATION COMPLETED'  ; call writeline
write(wline,'(a)')		; call writeline


return
end

subroutine read_MMp_matrix
!this subroutine reads the correlation matrix.  Data is stored in the data_module
use data_module
use write_module
use input_module
implicit none
integer(4)							:: i_ind
integer(4)							:: n_grp,i_grp

write(wline,'(a)')		; call writeline
write(wline,'(a)')'READ THE MMp MATRIX AND INDIVIDUAL INFORMATION'  ; call writeline
open(1,file=MMp_file,status='old',action='read')
read(1,*)n_ind

write(wline,'(a,i10)')'number of individuals                        :: ',n_ind	; call writeline
if(max_sig_ev == -1)then
	max_sig_ev=n_ind
endif
!n_ind=125
read(1,*)n_snp
if(n_snp_use == -1)then
	n_snp_use=n_snp
endif
write(wline,'(a,i10)')'number of SNPs on which this matrix is based :: ',n_snp	; call writeline
allocate(ind(n_ind)); ind%id='                    '
allocate(MMp(n_ind,n_ind))

do i_ind=1,n_ind
	read(1,*)ind(i_ind)%id(1:l_id),ind(i_ind)%sex,ind(i_ind)%dx,ind(i_ind)%grp,MMp(i_ind,:)
enddo

!print *,ind(n_ind)%id
!pause


close(1)
write(wline,'(a)'); call writeline
write(wline,'(a)')'COMPLETE Correlation (MMp) matrix read'	; call writeline
write(wline,'(a,i6)')'number of    males     :: ',count(ind%sex == 1)	; call writeline
write(wline,'(a,i6)')'number of  females     :: ',count(ind%sex == 2)	; call writeline
write(wline,'(a,i6)')'number of    cases     :: ',count(ind%dx == 2)	; call writeline
write(wline,'(a,i6)')'number of controls     :: ',count(ind%dx == 1)	; call writeline
n_grp =maxval(ind%grp)
do i_grp=1,n_grp
	write(wline,'(a,i3,a,i6)')'number in group :: ',i_grp,' :: ',count(ind%grp == i_grp)	; call writeline
enddo
write(wline,'(a)')	; call writeline
write(wline,'(a,120i6)')      '                      group ::        ',(i_grp,i_grp=1,n_grp)	; call writeline
write(wline,'(a,i6,1x,120i6)')'number of    male cases     :: ',count(ind%sex == 1 .and. ind%dx == 2), &
								(count(ind%sex == 1 .and. ind%dx == 2 .and. ind%grp == i_grp),i_grp=1,n_grp)	; call writeline
write(wline,'(a,i6,1x,120i6)')'number of  female cases     :: ',count(ind%sex == 2 .and. ind%dx == 2), &
								(count(ind%sex == 2 .and. ind%dx == 2 .and. ind%grp == i_grp),i_grp=1,n_grp)	; call writeline
write(wline,'(a,i6,1x,120i6)')'number of    male controls  :: ',count(ind%sex == 1 .and. ind%dx == 1), &
								(count(ind%sex == 1 .and. ind%dx == 1 .and. ind%grp == i_grp),i_grp=1,n_grp)	; call writeline
write(wline,'(a,i6,1x,120i6)')'number of  female controls  :: ',count(ind%sex == 2 .and. ind%dx == 1), &
								(count(ind%sex == 2 .and. ind%dx == 1 .and. ind%grp == i_grp),i_grp=1,n_grp)	; call writeline

return
end



subroutine read_input_data(input_file)
!subroutine reads the input parameters from the input file.  This information is stored in the input module
use input_module
use write_module
use parameter_module
implicit none
character(100)						:: input_file
integer(4)							:: i

!write(6,'(a)',advance='no')'give filename of file with input parameters :: '
!read(5,'(a)')input_file

open(10,file=input_file,status='old',action='read')
read(10,*)ext				!extension to use on the output files
read(10,*)data_dir			!directory where the data is stored, all files will be written to this directory
do i=1,150
	if(data_dir(i:i+1)== '/ ' .or. data_dir(i:i+1) == '\ ')then
		l_dir=i
		exit
	endif
enddo 

read(10,*)MMp_file				!Correlation matrix
read(10,*)exclusion_file		!Individuals to be excluded from this run
read(10,*)l_id					!length of the individuals identifier
read(10,*)min_cluster			!minimum cluster size
read(10,*)output				!amount of output generated 0=limited, 2=lots
read(10,*)max_sig_ev			!maximum number of eigenvalues to use (-1 = use based on data)
read(10,*)n_snp_use				!n_snp to use in calculation (-1 = used based on supplied data)
close(10)

open(16,file=data_dir(1:l_dir)//'GEM_log_'//ext//'.txt')

write(wline,'(a)')'==================================================='	; call writeline
write(wline,'(a)')'====           SPECTRAL GEM                    ===='	; call writeline
write(wline,'(a)')'==================================================='	; call writeline
write(wline,'(a)')											;call writeline
write(wline,'(a)')											;call writeline
write(wline,'(a,a)')'extension used on filenames :: ',ext	;call writeline
write(wline,'(a)')											;call writeline
write(wline,'(a)')											;call writeline


write(wline,'(a)')'input parameters'
write(wline,'(a)')      'name and location of the MMprime matrix :: '	;call writeline
write(wline,'(10x,a)')  MMp_file						;call writeline
write(wline,'(a)')      'name and location of the exclusion file :: '	;call writeline
write(wline,'(10x,a)')  exclusion_file				;call writeline
write(wline,'(a,i3,a,i10)')  'length of the identification field (max=',lid,')           :: ',l_id					;call writeline
write(wline,'(a)')																							;call writeline
write(wline,'(a,i10)')  'minimum cluster size                                   :: ',min_cluster			;call writeline
write(wline,'(a)')																							;call writeline
if(max_sig_ev == -1)then
	write(wline,'(a,i10)')  'Maximum number of significant dimension is based on the data and unlimited'	;call writeline
	write(wline,'(a)')																							;call writeline
else
	write(wline,'(a,i10)')  'Maximum number of significant dimensions to consider  :: ',max_sig_ev				;call writeline
	write(wline,'(a)')																							;call writeline
endif
if(n_snp_use == -1)then
	write(wline,'(a,i10)')  'SNP number used for calculations is based on the data    '							;call writeline
	write(wline,'(a)')																							;call writeline
else
	write(wline,'(a,i10)')  'SNP number used for calculations                     :: ',n_snp_use				;call writeline
	write(wline,'(a)')																							;call writeline
endif

return
end



subroutine write_minimum_distance_data
!subroutine that writes the data that is used for the minimum distance edit
use data_module
use input_module
use work_module
use write_module
implicit none
integer(4)						:: i_ind
integer(4)						:: m_ind

write(wline,'(a)')	; call writeline
write(wline,'(a)')'MINIMUM DISTANCES ARE STORED IN :: '	; call writeline
write(wline,'(5a)')'minimum_distance_'//ext//'.txt'	; call writeline

open(100,file=data_dir(1:l_dir)//'minimum_distance_'//ext//'.txt')
m_ind=count(ind%sub == 1)
do i_ind=1,m_ind
	write(100,'(a,3i2,2es13.5)')ind(list(i_ind))%id(1:l_id+1),ind(list(i_ind))%sex,ind(list(i_ind))%dx,ind(list(i_ind))%grp, &
							min_distance(list(i_ind),:)
enddo
close(100)

return
end

subroutine write_minimum_distance_R
!write the R program to use for deleting individuals based on minimum distances
use data_module
use input_module
use write_module
implicit none

write(wline,'(a)')	; call writeline
!write(wline,'(10a)')'R program to run after this step is :: ','case_control_distance_'//ext//'.R'	; call writeline

open(100,file=data_dir(1:l_dir)//'case_control_distance_'//ext//'.R')
write(100,'(a,a,a)')'data=read.table("minimum_distance_',ext,'.txt",header=FALSE)'
write(100,'(a,a,a)')'source("case_control_distance.R")'
write(100,'(a,a,a,a,a,a)')'write.table(c(exclude,"',ext,'"),file="exclude_dist_',ext, &
			'.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)'
close(100)

return
end


subroutine write_eigenvalue_decomposition(n_sig)
!this subroutine writes the information about the eigenvalue decomposition for this step
use data_module
use work_module
use input_module
use write_module
implicit none
integer(4)						:: n_sig
integer(4)						:: i_sig
integer(4)						:: i_ind
integer(4)						:: m_ind
integer(4)						:: dim_eval

write(wline,'(a)')	; call writeline
write(wline,'(a)')'SIGNIFICANT PART OF THE EIGENVALUE DECOMPOSITION FOR THIS STEP IS STORED IN :: '	; call writeline
write(wline,'(5a)')'significant_eigenvector_'//ext//'.txt'	; call writeline

dim_eval=size(eval)

!write the eigenvalues and vectors to a series of files
open(101,file=data_dir(1:l_dir)//'significant_eigenvalues_'//ext//'.txt')
do i_sig=1,n_sig
	write(101,'(i5,es13.5)')i_sig-1,eval(i_sig)
enddo
close(101)
open(100,file=data_dir(1:l_dir)//'significant_eigenvector_'//ext//'.txt')
m_ind=count(ind%sub == 1)
do i_ind=1,m_ind
	write(100,'(a,3i2,10000es13.5)')ind(list(i_ind))%id(1:l_id+1),ind(list(i_ind))%sex,ind(list(i_ind))%dx,ind(list(i_ind))%grp, &
							(evec(i_ind,i_sig),i_sig=1,n_sig)
!							(evec(i_ind,i_sig),i_sig=1,dim_eval)   !,dim_eval-n_sig+1,-1)
enddo
close(100)

return
end

subroutine write_matching_files(n_ev)
!this subroutine writes the file that is needed in the matching step
use write_module
use data_module
use input_module
use work_module
implicit none
integer(4)						:: n_ev
integer(4)						:: i_ind, m_ind
integer(4)						:: n_case, i_case
integer(4)						:: n_ctrl, i_ctrl
integer(4), dimension(:), allocatable	:: case_list
integer(4), dimension(:), allocatable	:: ctrl_list



write(wline,'(a)')	; call writeline
write(wline,'(a)')'WRITE THE FILE THAT IS USED FOR THE MATCHING STEP'	; call writeline

write(wline,'(a)')	; call writeline
write(wline,'(a,i3,a)')'calculate distances based on :: ',n_ev,' eigenvectors'; call writeline
call calculate_distance_NJW(n_ev,2,1)

n_case=count(ind%dx == 2 .and. ind%incl == 1)
n_ctrl=count(ind%dx == 1 .and. ind%incl == 1)
write(wline,'(a,i6)')'number of cases     :: ',n_case	; call writeline
write(wline,'(a,i6)')'number of controls  :: ',n_ctrl	; call writeline 
!make a list of cases and controls
allocate(case_list(n_case))
allocate(ctrl_list(n_ctrl))
n_case=0; n_ctrl=0
m_ind=count(ind%sub == 1)
do i_ind=1,m_ind
	if(ind(list(i_ind))%incl == 1)then
		if(ind(list(i_ind))%dx == 2)then
			n_case=n_case+1
			case_list(n_case)=i_ind
		else
			n_ctrl=n_ctrl+1
			ctrl_list(n_ctrl)=i_ind
		endif

	endif
enddo

write(wline,'(a)')	; call writeline
write(wline,'(a)')'Distance matrix is stored in :: '	; call writeline
write(wline,'(5a)')'distance_'//ext//'.matrix'	; call writeline
open(100,file=data_dir(1:l_dir)//'distance_'//ext//'.matrix')
if(n_ctrl >= n_case)then
	write(100,'(10000a)')(id_list(ctrl_list(i_ctrl))(1:l_id+1),i_ctrl=1,n_ctrl)
	do i_case=1,n_case
		write(100,'(a,100000es13.5)')id_list(case_list(i_case))(1:l_id+1),(dist(case_list(i_case),ctrl_list(i_ctrl)),i_ctrl=1,n_ctrl)
	enddo
else
	write(wline,'(a)')'number of cases greater than number of controls, the matrix is reversed'	; call writeline
	write(100,'(10000a)')(id_list(case_list(i_case))(1:l_id+1),i_case=1,n_case)
	do i_ctrl=1,n_ctrl
		write(100,'(a,10000es13.5)')id_list(ctrl_list(i_ctrl))(1:l_id+1),(dist(ctrl_list(i_ctrl),case_list(i_case)),i_case=1,n_case)
	enddo
endif

close(100)

deallocate(case_list)
deallocate(ctrl_list)



return
end

subroutine write_matching_R
!writes the R programs to use for FULL and PAIR matching
use write_module
use input_module
implicit none

write(wline,'(a)')	; call writeline
!write(wline,'(a)')'R PROGRAMS TO USE FOR MATCHING'	; call writeline

!write(wline,'(a)')	; call writeline
!write(wline,'(5a)')'Program to use for FULL matching :: ','full_matching_'//ext//'.R'	; call writeline

open(100,file=data_dir(1:l_dir)//'full_matching_'//ext//'.R')
write(100,'(10a)')'data=read.table("distance_',ext,'.matrix")'
write(100,'(10a)')'source("full_matched_set.R")'
write(100,'(10a)')'write.table(matches,file="fullmatch_',ext,'.txt",col.names=FALSE,quote=FALSE,sep=" ")'
close(100)


return
end


subroutine writeline
!this subroutine write a line of output to both the screen and the log-file.  It uses 
!the write module to transfer the data
use write_module
implicit none
integer(4)			:: i,j

do i=lline,1,-1
	if(wline(i:i) /= ' ')exit
enddo
i=i+1
!write(6,'(a)')wline(1:i)
write(16,'(a)')wline(1:i)

do j=1,i
	wline(i:i)=' '
enddo

return
end


