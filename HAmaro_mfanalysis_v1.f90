       program WAVE
	
       !---------------------------------------------------------------------------------------------------------------------------------------------
       ! Variables                                                                                                                                  !
       !---------------------------------------------------------------------------------------------------------------------------------------------
       
       implicit none

       !array allocations
       real(8), allocatable, dimension(:) :: signal
       real(8), allocatable, dimension(:) :: detail_coeff, approx_coeff

       !useful parameters, variables, indexes
       integer :: nnf
       integer(8) :: N1,N2,N3
       integer(8) :: nformat,ndatatype
       integer(8) :: nmax, L
       integer(8) :: NS, i, j, k, q, r ! -> indexes and length parameters
       integer, parameter :: file_id=13

       !variables related to written file names
       character(len=8) :: fmt, xi
       character (len = 16) :: file_name_16
       character (len = 13) :: structure_vec, u_vec, v_vec
       character (len = 100) :: ALIST
       character (len=13), allocatable, dimension(:) :: detail, approx, wavelet_leader_vec

       !variables to do with wavelet leaders determination
       integer :: length_j
       real(8) :: pos1, pos2, pos3

       !detail coefficients matrix, leader matrix and auxiliary matrix which is a portion of the detail coefficients matrix
       real(8), allocatable, dimension(:,:) :: djk_aux
       real(8), allocatable, dimension(:) ::  max_vec_top, pseudo_leaders, cjk_norm, djk_leader, djk_leader_q
       integer:: n_leader

       !variables to do with the scaling function
       integer :: qvalue, length_q !Scaling function parameter
       real(8) :: aux, delta_q, sum_result
       real(8), allocatable, dimension(:,:) ::  U_jq, V_jq, S_jq   !scaling function matrix
       real(8), allocatable, dimension(:) :: d_spec, h_spec, X, wj, zeta_q
       integer, allocatable, dimension(:):: leader_count, q_vec
       real(8) :: bj
       integer :: j1,j2, j1_input, j2_input

       !cumulants
       real(8), allocatable, dimension(:,:) :: CP
       real(8) :: v0, v1, v2, m1, m2, m3
       real(8), dimension(3) :: cumul
       ! writing files variables
       integer :: p_coeffs, p_leaders, p_sgram, p_results, p_struct, p_uv

       ! clock
       integer:: t1, t2
       integer :: clock_rate
       integer :: clock_max

       !!!!!!!! OPENMP !!!!!!!!!!!!!!!!
       integer :: max_num_threads, nps
       integer :: omp_get_max_threads
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !--------------------------------------------------------------------------------------------------------------------
       nnf=100

       !file_name='HWData_00000.out'
       !file_name='sinedd_00000.out'
       !file_name='cosine_00000.out'
       !file_name='rlix_nv_0002.scl'

       call system_clock ( t1, clock_rate, clock_max )
       
       !---------------------------------------------------------------------------------------------------------------------------------------------
       ! Reading input parameters from file 'input_multifractal.tex'                                                                                !
       !---------------------------------------------------------------------------------------------------------------------------------------------
       
       write(*,*)
       write(*,*) 'Reading parameters'
       write(*,*)
       open(nnf,file='input_mfanalysis_v1.tex',status='unknown')

       read(nnf,*) file_name_16
       write(*,*) 'File name: ',file_name_16
       read(nnf,*) nps
       write(*,*) 'Number threads: ',nps
       
       read(nnf,*) nformat
       
       if (nformat.eq.1) then
               write(*,*) 'Data format: ASCII'
       else if (nformat.eq.2) then
               write(*,*) 'Data format: Binary'
       endif        
       
       read(nnf,*) ndatatype
       if (ndatatype.eq.1) then
              write(*,*) 'Data Type: 1D Vector'
       else if (ndatatype.eq.3) then
              write(*,*) 'Data Type: 3D Box [N1xN2xN3]'
       endif 
       
       read(nnf,*) ALIST
       
       print*,''
       
       print*, 'Dimensions'
       read(nnf,*) N1
       write(*,*) 'N1 ',N1
       read(nnf,*) N2
       write(*,*) 'N2: ',N2
       read(nnf,*) N3
       write(*,*) 'N3: ',N3
       write(*,*)
       read(nnf,*) ALIST
       read(nnf,*) p_coeffs
       read(nnf,*) p_leaders
       read(nnf,*) p_sgram
       read(nnf,*) p_results
       read(nnf,*) p_struct
       read(nnf,*) p_uv
       read(nnf,*) ALIST
       read(nnf,*) j1_input
       read(nnf,*) j2_input

       close(nnf)
       print*,'Quantities to be printed:'
       print*,'coefficients|','    leaders|','   scalogram|','    results|','   struct func|', '   U and V|'
       print*,p_coeffs,p_leaders,p_sgram,p_results, p_struct, p_uv
       call OMP_SET_NUM_THreadS(nps)

       max_num_threads = omp_get_max_threads()
       
       print*,''
       if(nps>max_num_threads)then
              print*, 'WARNING --- SELECTED N PROCESSORS =',nps,'MAX PROC AVAIL =',max_num_threads
              nps = max_num_threads
       end if

       !---------------------------------------------------------------------------------------------------------------------------------------------
       ! Reading data from file specified in input file 'input_multifractal.tex'                                                                    !
       !---------------------------------------------------------------------------------------------------------------------------------------------
       
       call read_array_all(file_name_16,nformat,ndatatype, N1,N2,N3, signal,nmax)
       write(*,*)

       !obtaining the moments vector-----------------------------------------------------------------------------------------------------------
       
       qvalue=5 !--------------------------> inputs
       delta_q=1!--------------------------> inputs
                
       !length_q=int(qvalue/delta_q) +1
       length_q=int(2*(qvalue/delta_q)) +1 !-> if considered negative moments q < 0
                
       !moments vector
       q_vec= (/ (-qvalue + (i-1)*delta_q, i=1,length_q)/) !-> if considered negative moments < 0
       
       !----------------------------------------------------------------------------------------------------------------------------------------

       !signal length
       L=nmax
       !number of scales
       NS=nint(log(real(L))/log(2D0))

       if (j1_input.eq.0) then 
               j1=3    !----------------------------------------------------> first scale considered
       else
               j1=j1_input
       end if

       if (j2_input.eq.0) then 
               j2=NS-3    !----------------------------------------------------> final scaled comsidered
       else
               j2=j2_input
       end if
       

       write(*,*) 'Number of used points, N =',nmax
       print*,'Number of Scales =',NS
       print*,'Relevant scale boundaries: '
       print*, 'J1 =',j1
       print*, 'J2 =',j2
       
       !create array with name of the dwt resultant files--------------------------------------------------------------------
       allocate(detail(NS))
       allocate(approx(NS))
       allocate(wavelet_leader_vec(NS))

       !create the output files' name arrays
       fmt = '(I3.3)' ! an integer of width 3 with zeros at the left

       do i=1,NS
              write (xi,fmt) i ! converting integer to string using a 'internal file'
              detail(i)='dwt_detail'//trim(xi)
              approx(i)='dwt_approx'//trim(xi)
              wavelet_leader_vec(i)='wav_leader'//trim(xi)
       end do
      

       ! open file for scalogram
       if(p_sgram.eq.1) then
               open(101,file='scalogram.plt',status='unknown')
               write(101,*) 'VARIABLES = "Time", "Frequency", "Map"'
               write(101,*) 'ZONE I=',2**(NS-1),'J=',NS-1,'F=POINT'
       endif
       
       ! Cascade
       allocate(cjk_norm(L/2))
       allocate(djk_leader(L/2))
       allocate(pseudo_leaders(L/2))
       allocate(max_vec_top(L/2))
       allocate(leader_count(j2+1))
       allocate(U_jq(NS-1,length_q))
       allocate(V_jq(NS-1,length_q))
       allocate(djk_leader_q(L/2))
       allocate(S_jq(NS-1,length_q))
       allocate(wj(NS-1))
       allocate(CP(3,NS-3))
       
       cjk_norm=0.
       aux=10e8
       sum_result=0.
       djk_leader=0.
       leader_count=0
       wj=0.
       CP=0.
       
       print*,''
       print*, 'Calculating Wavelet leaders'
       
       !------------------------------------
       allocate(X(NS-1))
       ! X is the scales vector       
       X=0.        
       X = (/ ( i, i=1,NS-1)/)     
       
! MAIN CYCLE ----------------------------------------------------------------------------------------------------------------------------------       
       do j=1,NS-1
              
              length_j=L/(2**j)


              !--------------------------------------------------------------------------------------------------------------------
              !  Performing the discrete wavelet transform to the signal                                                          !
              !--------------------------------------------------------------------------------------------------------------------  
              
              call dwt(signal,approx_coeff,detail_coeff) !dwt execute
              signal=approx_coeff
              
              !build the normalized coefficient matrix cjk_norm
              do k=1,length_j-2
              cjk_norm(k+1)=2**(-j*0.5)*abs(detail_coeff(k))
              end do

              !eliminate border effect by making boarders inf
              cjk_norm(1)=10e10
              cjk_norm(length_j)=10e10

              !--------------------------------------------------------------------------------------------------------------------
              !  Obtained the wavelet leaders                                                                                     !
              !--------------------------------------------------------------------------------------------------------------------              
              max_vec_top=0.
              pos1=0.
              pos2=0.
              pos3=0.
              
              if (j.eq.1) then
              
                     allocate(djk_aux(3,L/2))
                     djk_aux=0.
    
                     djk_aux(1,:)=cjk_norm(1:length_j)
                     djk_aux(2,1:length_j-1)=cjk_norm(2:length_j)
                     djk_aux(3,1:length_j-2)=cjk_norm(3:length_j)
                     pseudo_leaders(:)=cjk_norm(:)
    
                     djk_leader=maxval(djk_aux,DIM=1)
    
                     deallocate(djk_aux)
    	
              elseif (j.ge.2) then
                     !djk_aux is a auxiliary matrix that has previous scale leaders
                     
                     allocate(djk_aux(3,length_j))
                     djk_aux=0.
    
                     djk_aux(1,:)=cjk_norm(1:length_j)
                     djk_aux(2,:)=pseudo_leaders(1:length_j*2:2)
                     djk_aux(3,:)=pseudo_leaders(2:length_j*2:2)
    
                     pseudo_leaders=maxval(djk_aux,DIM=1)
                   
    
                     do r=1,length_j-2
                            if (r.eq.(length_j-1)) then
                                   pos1=pseudo_leaders(r)
                                   pos2=pseudo_leaders(r+1)
                                   djk_leader(r)=MAX(pos1,pos2)
                           
                                   elseif (r.eq.(length_j)) then
                                   pos1=pseudo_leaders(r)
                                   djk_leader(r)=pos1
                           
                                   else
                                   pos1=pseudo_leaders(r)
                                   pos2=pseudo_leaders(r+1)
                                   pos3=pseudo_leaders(r+2)
                                   djk_leader(r)=MAX(pos1,pos2,pos3)
                            endif
                     end do

                     do k=1,length_j
                            if (djk_leader(k)>aux) djk_leader(k)=0
                     end do
              
                     deallocate(djk_aux)
               endif

               !clean the entries where infinity lies
               do k=1,length_j
                      if (djk_leader(k)>aux) djk_leader(k)=0
               end do

              djk_leader=pack(djk_leader(:),djk_leader(:)>0)
              
              n_leader = count(djk_leader(:)/=0) !counting the leaders in each scale
              
              if (n_leader.eq.0) exit !exit clause when leaders is empty
            
              !write leaders into files
              if ((p_leaders.eq.1).and.(n_leader.ne.0)) then
                     open(101, file=wavelet_leader_vec(j), action="write")
                     call writex(101,djk_leader(1:n_leader))
                     close(101)
              endif
              
              !--------------------------------------------------------------------------------------------------------------------
              ! Obtaining U and V fields and structure function S                                                                 !
              !--------------------------------------------------------------------------------------------------------------------

              
              !$OMP PARALLEL &
              !$OMP DEFAULT(none) &
              !$OMP PRIVATE(q,djk_leader_q,sum_result) &
              !$OMP SHARED(djk_leader,q_vec,V_jq,U_jq,S_jq,length_q,n_leader,j2,j)
            
              if (j.le.j2) then
              !$OMP DO &
              !$OMP SCHEDULE(static)
                     do q=1,length_q

                            djk_leader_q=djk_leader(1:n_leader)**q_vec(q) !leader matrix to the power of the scaling moment q
			
	                    sum_result=sum(djk_leader_q(1:n_leader),1) !sum of the leaders for each scale j
                            
                            !$OMP CRITICAL
                            V_jq(j,q)=(1/sum_result)* &
	                    sum(djk_leader_q(1:n_leader)*(log(djk_leader(1:n_leader))/log(2.)),1)
			
		             U_jq(j,q)=(log(real(n_leader))/log(2.))+ &
                             sum((djk_leader_q(1:n_leader)/sum_result)* &
                             (log(djk_leader_q(1:n_leader)/sum_result)/log(2.)),1)
                            
                            S_jq(j,q)=log(sum(djk_leader_q(1:n_leader),1)/n_leader)/log(2.)
                            !$OMP END CRITICAL
                      enddo
               !$OMP END DO
               
               endif
               !$OMP END PARALLEL
               !-------------------------------------------------------------------------------------------------------------------
               ! Calculating weights for linear regressions                                                                       !
               !-------------------------------------------------------------------------------------------------------------------
                      bj=n_leader
              
                      !bj=1
                 
                      v0=sum((X(j1:j2)**0)*bj,1)
                      v1=sum((X(j1:j2)**1)*bj,1)
                      v2=sum((X(j1:j2)**2)*bj,1)
                
                      wj(j)= bj*((v0*j-v1)/(v0*v2-(v1**2)))	!weight vector for linear regression

               !-------------------------------------------------------------------------------------------------------------------
               ! Calculating the cumulants                                                                                        !
               !-------------------------------------------------------------------------------------------------------------------
               
                      m1=sum(log(djk_leader(1:n_leader)),1)/n_leader
                      CP(1,j)=m1
                      m2=(sum(log(djk_leader(1:n_leader))**2,1)/n_leader)
                      CP(2,j)=m2-(m1**2)
                      m3=(sum(log(djk_leader(1:n_leader))**3,1)/n_leader)
                      CP(3,j)=m3-3*(m2*m1)+2*m1**3

               !-------------------------------------------------------------------------------------------------------------------
               ! Printing coefficients and scalogram                                                                              !
               !-------------------------------------------------------------------------------------------------------------------
               if (p_coeffs.eq.1) then
            
                      !write approximation coefficients into files
                      open(file_id, file=approx(j), action="write")
                      call writex(file_id, approx_coeff)
                      close(file_id)
            
                      !write detail coefficients into files
                      open(file_id, file=detail(j), action="write")
                      call writex(file_id, cjk_norm(:))
                      close(file_id)
               endif
            
               !write scalogram for tecplot
               if (p_sgram.eq.1) then
            
                      do k=1,L/(2**j)
                             if (cjk_norm(k).ge.aux) then
                                    write(101,*) (2**(j-1)*(k-1)+1),j,0D2
                             else
                                    write(101,*) (2**(j-1)*(k-1)+1),j,cjk_norm(k)
                             endif
                      enddo
               endif
       
       djk_leader=0.
       
       end do
       
!-----------------------------------------------------------------------------------------------------------------------------------------------------
                    
       close(101)
       deallocate(detail)
       deallocate(approx)
       deallocate(signal)
       deallocate(cjk_norm)
       deallocate(wavelet_leader_vec)
       deallocate(max_vec_top)
       deallocate(pseudo_leaders)
       print*,'(Done)'         
                
       !-------------------------------------------------------------------------------------------------------------------
       ! Obtaining the spectrum (h and d) the scaling exponenent and cumulants using linear regression                    !
       !-------------------------------------------------------------------------------------------------------------------                

                
       allocate(d_spec(length_q))
       allocate(h_spec(length_q))
       allocate(zeta_q(length_q))
       
       do q=1,length_q
       
              !Chhabra's method
              h_spec(q)=sum(wj(j1:j2)*V_jq(j1:j2,q))
              d_spec(q)=1+sum(wj(j1:j2)*U_jq(j1:j2,q))
              !---------------------------------------------
              zeta_q(q)=sum(wj(j1:j2)*S_jq(j1:j2,q))
       end do     
       
       CP=CP/log(2.)
       cumul(1)=sum(CP(1,j1:j2)*wj(j1:j2),1)
       cumul(2)=sum(CP(2,j1:j2)*wj(j1:j2),1)
       cumul(3)=sum(CP(3,j1:j2)*wj(j1:j2),1)
      
    
       !------------------------------------------------------------------------------------------------------------------
       ! Printing quantities on screen                                                                                   !
       !------------------------------------------------------------------------------------------------------------------
       print*,''
       print*,'Results:'
       print*,'------------------------------------------------------'
       print*,'d_spec ----->', d_spec
       print*,'------------------------------------------------------'
       print*,'h_spec ----->', h_spec
       print*,'------------------------------------------------------'
       print*,'zeta_q ----->', zeta_q
       print*,'------------------------------------------------------'
       print*,'cumult ----->', cumul
       print*,'------------------------------------------------------'
    
       !------------------------------------------------------------------------------------------------------------------
       ! Writing quantities into files                                                                                   !
       !------------------------------------------------------------------------------------------------------------------
       
       if(p_results.eq.1) then
       !write d_spec into file
       open(file_id, file="d_spec.txt", action="write")
       call writex(file_id, d_spec)
       print*,'Written file, "d_spec.txt" '
       close(file_id)
    
       !write h_spec into file
       open(file_id, file="h_spec.txt", action="write")
       call writex(file_id, h_spec)
       print*,'Written file, "h_spec.txt" '
       close(file_id)
    
       !write scaling coefficients zeta(q) into file
       open(file_id, file="zeta_q.txt", action="write")
       call writex(file_id, zeta_q)
       print*,'Written file, "zeta_q.txt" '
       close(file_id)
    
       !write cumulants into file
       open(file_id, file="cumul.txt", action="write")
       call writex(file_id, cumul)
       print*,'Written file, "cumul.txt" '
       close(file_id)
       endif
       
       !print structure function into files 
       if(p_struct.eq.1) then
              do q=1,length_q
                     if (q_vec(q).lt.0) then
                            fmt = '(I1.0)'
                            write (xi,fmt) abs(q_vec(q))   
                            structure_vec='structf_q_0'//trim(xi)
                     elseif (q_vec(q).gt.0) then
                            fmt = '(I1.0)'
                            write (xi,fmt) q_vec(q)
                            structure_vec='structf_q_'//trim(xi)
                     endif
                     open(file_id, file=structure_vec, action="write")
                     call writex(file_id, pack(S_jq(:,q),S_jq(:,q)/=0))
              enddo
              print*,'Written files, "structf_j" '
              close(file_id)
       endif
       
       !print U,V fields into files
       if(p_uv.eq.1) then
              do q=1,length_q
                     if (q_vec(q).lt.0) then
                            fmt = '(I1.0)'
                            write (xi,fmt) abs(q_vec(q))   
                            u_vec='U_q_0'//trim(xi)
                            v_vec='V_q_0'//trim(xi)
                     elseif (q_vec(q).gt.0) then
                            fmt = '(I1.0)'
                            write (xi,fmt) q_vec(q)
                            u_vec='U_q_'//trim(xi)
                            v_vec='V_q_'//trim(xi)
                     else
                            u_vec='U_q_00'
                            v_vec='V_q_00'
                            
                     endif
                     open(file_id, file=u_vec, action="write")
                     call writex(file_id, pack(U_jq(:,q),U_jq(:,q)/=0))
                     open(file_id, file=v_vec, action="write")
                     call writex(file_id, pack(V_jq(:,q),V_jq(:,q)/=0))
              enddo
              print*,'Written files, "structf_j" '
              close(file_id)
       endif
       !---------------------------------------------------------------------------------------------------------------------------------------------    

       deallocate(d_spec)
       deallocate(h_spec)
       deallocate(U_jq)
       deallocate(V_jq)
       deallocate(wj)
       deallocate(leader_count)
       deallocate(X)
       deallocate(S_jq)
       deallocate(djk_leader_q)
       deallocate(q_vec)
       
       
       
       call system_clock ( t2, clock_rate, clock_max )
       print*,''
       write ( *, * ) 'Elapsed real time (s) = ', real ( t2 - t1 ) / real ( clock_rate )
      
!---------------------------------------------------------------------------------------------------------------------------------------------------!      
       
       !--------------------------------------------------------------------------------------------------------------------------------------------!
       !---------------------------------------------------------------------------------------------------------------------------------------------
       ! Subroutines                                                                                                                                !
       !--------------------------------------------------------------------------------------------------------------------------------------------- 
       !--------------------------------------------------------------------------------------------------------------------------------------------!
       
       contains
       
       !---------------------------------------------------------------------------------------------------------------------------------------------
       subroutine dwt(input_array,trans1,trans2)
    
       implicit none
    
       integer::L
       real(8),allocatable,dimension(:):: temp1, temp2
       real(8),allocatable,dimension(:),intent(inout):: input_array
       real(8),allocatable,dimension(:),intent(out):: trans1, trans2
    
       !scaling function coefficients for Daubechies 4 aka approx coefficients
       real(8),parameter::C0=0.4829629131445341, C1=0.8365163037378079, C2=0.2241438680420134, C3=-0.1294095225512604
    
       !wavelet function coefficients for Daubechies 4 aka detail coefficients
       real(8),parameter::G0=C3, G1=-C2, G2=C1, G3=-C0
    
       L=size(input_array)
    
       if (L>1) then
              allocate(temp1(L/2))
              allocate(temp2(L/2))
              allocate(trans1(L/2))
              allocate(trans2(L/2))
        
              !approximation coefficients
              temp1(1:L/2-1)=C0*input_array(1:L-3:2)+C1*input_array(2:L-2:2)+C2*input_array(3:L-1:2)+C3*input_array(4:L:2)
              temp1(L/2)=C0*input_array(L-1)+C1*input_array(L)+C2*input_array(1)+C3*input_array(2)
        
              !detail coefficients
              temp2(1:L/2-1)=G0*input_array(1:L-3:2)+G1*input_array(2:L-2:2)+G2*input_array(3:L-1:2)+G3*input_array(4:L:2)
              temp2(L/2)=G0*input_array(L-1)+G1*input_array(L)+G2*input_array(1)+G3*input_array(2)
        
              trans1=temp1
              trans2=temp2
        
              deallocate(temp1)
              deallocate(temp2)
       endif
        
        
       end subroutine dwt
      !----------------------------------------------------------------------------------------------------------------------------------------------
       subroutine writex(file_id, vec)
       real(8), dimension(:), intent(in) :: vec
       integer, intent(in) :: file_id
       integer :: i
       write(file_id, "(F30.16)") ( vec(i), i=lbound(vec,1), ubound(vec,1))
    
       end subroutine writex

      !----------------------------------------------------------------------------------------------------------------------------------------------
       subroutine read_array_all(file_name_16,nformat,ndatatype, N1,N2,N3, signal,nmax)
       implicit none
	
       !variables
       integer(8) :: nnf
       integer(8) :: nmax,l
       integer(8) :: i,j,k
       integer(8) :: N1,N2,N3
       integer(8) :: nformat,ndatatype


       integer(8) :: filesize,reclength,recidx

       real(8), allocatable, dimension(:) :: signal
       real(4), allocatable, dimension(:,:,:) :: aux3D

       character (len = 16) :: file_name_16

       nnf=100
    
       allocate(aux3D(N1,N2,N3))
    
       ! compute nmax
       if (ndatatype.eq.1) then
              nmax=N1
       elseif (ndatatype.eq.2) then
              nmax=N1*N2
       elseif (ndatatype.eq.3) then
              nmax=N1*N2*N3
       endif
       
       !write(*,*) 'nmax = ',nmax
        
       allocate(signal(nmax))
        
       ! reading raw data into signal
       if (ndatatype.eq.1) then
            
              write(*,*) 'Reading file ',file_name_16
              if (nformat.eq.1) then
                     open(unit=nnf,file=file_name_16,status='unknown')
                     do i=1,nmax
                            l=i
                            read(nnf,*) signal(l)
                     enddo
                     close(nnf)
              endif
              
              write(*,*) '(done)'
                
       elseif (ndatatype.eq.3) then
       
       
              write(*,*) 'Reading 3D file ',file_name_16
              inquire (file = file_name_16, size = filesize)
              inquire (iolength=reclength)            &
                     ((aux3D(i,j,1), i=1,N1), j=1,N2)
              open (unit=nnf,FILE=file_name_16,STATUS='OLD',form='UNFORMATTED',access="stream")
              
              do k=1,N3
                      recidx=(k-1)*reclength+1
                      read(unit=nnf,POS=recidx)           &
                      ((aux3D(i,j,k), i=1,N1), j=1,N2)
             enddo
             close(nnf)
             write(*,*) '(done)'
              
             
             l=0
              
              do k=1,N3
                     do j=1,N2
                            do i=1,N1
                                    l=l+1
                                    signal(l)=aux3D(i,j,k)
                            enddo
                     enddo
              enddo
       endif
                
                
       deallocate(aux3D)
                
       end subroutine read_array_all
      !---------------------------------------------------------------------------------------------------------------------------------------------
end program WAVE
