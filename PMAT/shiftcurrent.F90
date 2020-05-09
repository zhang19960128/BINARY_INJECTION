program shift
  USE constants
  USE flinal_tools, ONLY : ridentity, iidentity
  USE mp_global, ONLY : nproc, proc_id, mp_start, mp_end, mp_barrier, mp_comm_all, mp_bcast
  USE read_tools, ONLY : header_kl, read_wfc_kl
  USE io_wfc, ONLY : trans_qewfc
  USE op_tools, ONLY : vector_op, kpt_ops, lattice, op_blocks_kl, release_matrices_data, dump_op_kl, p_block,vector_op_two,s_block,scalar_op
  USE es_tools, ONLY : check_kpt_id_list, kpt_init, es_header, release_kpt_kl, inner_prod
  USE sc_tools, ONLY : shift_current_data, shift_kpt, dump_sg_kl, dump_sc_kl

  USE sc_tools, ONLY : energy_grid_ev, release_sc_data, &
      dielectric_spectrum_kl, collect_spectrum_di, write_di_spectrum_kl, &
      shift_current_spectrum2_kl, collect_spectrum_sc, collect_spectrum_sc2, write_sc_spectrum_kl, &
      shift_vector_spectrum_kl, collect_spectrum_sv, write_sv_spectrum_kl, &
      calculate_glass_spectrum, calculate_glass_spectrum2, write_glass_spectrum_kl, &
      sg_spectrum_2omega_kl, collect_spectrum_sg, write_sg_spectrum_kl

  USE control, ONLY : dirname, kpts_type, data_format, cutoff, broadwidth, &
                      block_tol, sv_tol, sg_tol, resolution, sg_calc

  USE es_tools, ONLY : get_kpt_grid_kl, get_kpt_id_list, get_kpt_id_link
  USE sc_readin, ONLY : read_input
  USE io_global, ONLY : stdout, ionode
  USE version, ONLY : version_number, min_qe_version
  use MPI
implicit none
INTEGER :: ierr
type(es_header),pointer                                                     :: es_main_kl, es_test
type(kpt_ops)                                                               :: matrices
type(shift_current_data)                                                    :: sc
type(vector_op)                                                             :: p_matrix_block
type(scalar_op)                                                         :: inner_matrix_block

real(dp),dimension(:,:,:,:,:), allocatable                                  :: sc_spec
real(dp),dimension(:,:,:,:,:,:), allocatable                                :: sg_spec
real(dp),dimension(:,:,:,:,:), allocatable                                  :: pre_sc_spec
real(dp),dimension(:,:,:,:,:), allocatable                                  :: di_spec
real(dp),dimension(:,:,:), allocatable                                      :: sv_spec
double precision,dimension(:),allocatable                                :: kpointmomentum       
real(dp),dimension(:), allocatable                                          :: egrid
real(dp),dimension(3)                                                       :: kpoint_define
complex(dpc)                                                                :: overlap,momentum_k
complex(dpc),dimension(3)                                                   :: momentum,momentum_periodical
logical                                                                     :: store_wf=.true., test_mod=.false., nonmp

integer                                                                     :: ikptset,jkpt,iblock,jblock,ibnd,jbnd,ipw, k(3), kcoord(3), nspin, ispin, output_format, narg,iband1,iband2,kpointindex
integer                                                                     :: mpirank,mpisize,buffsize,mpifile,kpointid,tempinteger,totalbands,oneDindex
integer(kind=MPI_OFFSET_KIND)                                               :: offset               
integer,dimension(3)                                                        :: buffinteger
integer,dimension(10)                                                       :: buffreal      
integer                                                                     :: status(mpi_status_size)
double precision                                                            :: tempreal           
  !
  ! ... Initialize MPI
  CALL mp_start
  !
  IF (ionode) THEN
    WRITE(stdout,'(A,A)') "VERSION ", version_number
    WRITE(stdout,'(A,A)') "Require QE version >= ", min_qe_version
    WRITE(stdout,*)
  END IF
  !
  ! ... Read control parameters from input
  CALL read_input
  !
  ! ... decide output format for electronic structure calc
  CALL decide_input(output_format)
  !
  ! ... Read system details from QE output
  es_main_kl=>header_kl(output_format)
  !
  ! ... Done with initialization
  CALL mp_barrier(mp_comm_all)
  !
  ! ... Run main calculation
  call mpi_file_open(MPI_COMM_WORLD,"pmat.binary",MPI_MODE_WRONLY +MPI_MODE_CREATE,MPI_INFO_NULL,mpifile,ierr)
  DO ikptset=1,size(es_main_kl%kpt_id_link,3)
    kpointid=(ikptset-1)*size(es_main_kl%kpt_id_link,2)+proc_id
    DO ispin=1,es_main_kl%nspin
      !
!           build es_header
            call flush(6)
            call read_wfc_kl( es_main_kl,ikptset,ispin, output_format)
            call mp_barrier(mp_comm_all)
!           build op_block data
            matrices=op_blocks_kl(es_main_kl,cutoff,sg_calc)
  !             Do ibnd=1,size(es_main_kl%kpt(ikptset)%block(iblock)%wf)
  !              Do jbnd=1,size(es_main_kl%kpt(ikptset)%block(jblock)%wf)
  !              write (22,'(I6,6(E20.10))') ikptset, real(p_matrix_block%v(ibnd,jbnd,:)),aimag(p_matrix_block%v(ibnd,jbnd,:))
  !              end do
  !             end do
  !             Do ibnd=1,size(es_main_kl%kpt(0)%block(iblock)%wf)
  !              Do jbnd=1,size(es_main_kl%kpt(0)%block(jblock)%wf)
  !          !    write (22,'(I6,6(E20.10))') ikptset, real(matrices%block(iblock,jblock)%p%v(ibnd,jbnd,:)),aimag(matrices%block(iblock,jblock)%p%v(ibnd,jbnd,:))
  !              end do
  !             end do
             iband1=1
             iband2=1
             totalbands=es_main_kl%nband
             if(allocated(es_main_kl%kpt(kpointindex)%block)) then
             if(not(allocated(kpointmomentum))) then
               allocate(kpointmomentum(totalbands*totalbands*10))
             end if
             kpointindex=0
             kpoint_define=es_main_kl%kpt(kpointindex)%g%k
             Do iblock=1,es_main_kl%kpt(kpointindex)%nblock
                Do ibnd=1,size(es_main_kl%kpt(kpointindex)%block(iblock)%wf)
                   DO jblock=1,es_main_kl%kpt(kpointindex)%nblock
                     Do jbnd=1,size(es_main_kl%kpt(kpointindex)%block(jblock)%wf)
                      p_matrix_block=p_block(es_main_kl%kpt(kpointindex)%block(iblock)%wf,es_main_kl%kpt(kpointindex)%block(jblock)%wf,es_main_kl%lattice%gprimd)
                      inner_matrix_block=s_block(es_main_kl%kpt(kpointindex)%block(iblock)%wf,es_main_kl%kpt(kpointindex)%block(jblock)%wf)
                           momentum_periodical=p_matrix_block%v(ibnd,jbnd,:)
                           momentum_k=inner_matrix_block%s(ibnd,jbnd)
                           momentum=momentum_periodical
                           momentum(1)=momentum(1)+momentum_k*kpoint_define(1)
                           momentum(2)=momentum(2)+momentum_k*kpoint_define(2)
                           momentum(3)=momentum(3)+momentum_k*kpoint_define(3)
                           oneDindex=((iband1-1)*totalbands+(iband2-1))*10
                           Do tempinteger=1,3
                           kpointmomentum(oneDindex+tempinteger)=real(momentum(tempinteger))
                           kpointmomentum(oneDindex+tempinteger+3)=aimag(momentum(tempinteger))
                           end do
                           kpointmomentum(oneDindex+7)=es_main_kl%kpt(kpointindex)%block(iblock)%energy
                           kpointmomentum(oneDindex+8)=es_main_kl%kpt(kpointindex)%block(jblock)%energy
                           kpointmomentum(oneDindex+9)=es_main_kl%kpt(kpointindex)%block(iblock)%filling
                           kpointmomentum(oneDindex+10)=es_main_kl%kpt(kpointindex)%block(jblock)%filling
                           iband2=iband2+1
                           if (iband2==es_main_kl%kpt(kpointindex)%nband+1) then
                           iband1=iband1+1
                           iband1=mod(iband1-1,es_main_kl%kpt(kpointindex)%nband)+1
                           endif
                           iband2=mod(iband2-1,es_main_kl%kpt(kpointindex)%nband)+1
                     end do
                   end do
                  end do
               end do
            offset=kpointid*totalbands*totalbands*10*sizeof(tempreal)
            write(*,*) "the offset is: ",offset,"the total number of numbers I need to write is: ",sizeof(kpointmomentum)
            call MPI_File_Write_at(mpifile,offset,kpointmomentum,totalbands*totalbands*10,MPI_DOUBLE_PRECISION,status,ierr)
            !if doing a kpoint path, just print out and end here
            end if
            if (kpts_type=='path') then
              call dump_op_kl(matrices,int(es_main_kl%nspin/2.0)+ispin,cutoff)
              call release_kpt_kl(es_main_kl)
              cycle
            end if
            call release_kpt_kl(es_main_kl)
            call mp_barrier(mp_comm_all)
            CALL release_matrices_data(matrices)
            if (proc_id==0) write(6,*) 'deallocated sc_date for ikptset', ikptset, 'ispin=', ispin
        end do  !loop over ispin
        call mp_barrier(mp_comm_all)
    end do
  call mpi_file_close(mpifile,ierr)
  call mp_barrier(mp_comm_all)
  write(*,*) "finished without error"
  CALL mp_end
END PROGRAM
