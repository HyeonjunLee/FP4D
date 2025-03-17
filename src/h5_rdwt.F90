module H5_RDWT

USE HDF5 ! This module contains all necessary modules

! in this module, all k2 means T and all lambda means loglambda
CONTAINS

	! read velocity grid for test interpolation
 	! 
	subroutine h5_vread(vrMatrix, vrpMatrix, vzMatrix, vzpMatrix)  

		IMPLICIT NONE
		CHARACTER(LEN=17) 		:: filename = "UperpperpCheck.h5"
		CHARACTER(LEN=8)	     	:: vrname = "vrMatrix"    	
		CHARACTER(LEN=9)	     	:: vrpname = "vrpMatrix"    
		CHARACTER(LEN=8)	     	:: vzname = "vzMatrix"
		CHARACTER(LEN=9)	     	:: vzpname = "vzpMatrix"    			
		INTEGER		, PARAMETER	:: checkdim = 140
		INTEGER 			:: error 					! Error flag
		INTEGER(HID_T) 			:: file_id		 			! file handlers
		INTEGER(HID_T) 			:: vrid, vrpid, vzid, vzpid			! Nint data handles
		INTEGER(HSIZE_T), DIMENSION(1)  :: dims = (/checkdim/)				! check data size buffer
		REAL(KIND=8)    , DIMENSION(checkdim)	, intent(out) :: vrMatrix, vrpMatrix, vzMatrix, vzpMatrix	! read velocity data buffer

		! Initialize FORTRAN interface.
		! 
		CALL h5open_f(error)
		CALL h5fopen_f (filename, H5F_ACC_RDONLY_F, file_id, error)

		! Open an vr dataset and read
		! 
		CALL h5dopen_f(file_id, vrname, vrid, error)
		CALL h5dread_f(vrid, H5T_IEEE_F64LE, vrMatrix, dims, error)
		CALL h5dclose_f(vrid, error)

		! Open an vrp dataset and read
		! 
		CALL h5dopen_f(file_id, vrpname, vrpid, error)
		CALL h5dread_f(vrpid, H5T_IEEE_F64LE, vrpMatrix, dims, error)
		CALL h5dclose_f(vrpid, error)

		! Open an vz dataset and read
		! 
		CALL h5dopen_f(file_id, vzname, vzid, error)
		CALL h5dread_f(vzid, H5T_IEEE_F64LE, vzMatrix, dims, error)
		CALL h5dclose_f(vzid, error)

		! Open an vzp dataset and read
		! 
		CALL h5dopen_f(file_id, vzpname, vzpid, error)
		CALL h5dread_f(vzpid, H5T_IEEE_F64LE, vzpMatrix, dims, error)
		CALL h5dclose_f(vzpid, error)

		! Close the file.
		! 
		CALL h5fclose_f(file_id, error)
		CALL h5close_f(error)

	end subroutine h5_vread

	! 
	! read all table value for 3D data (rdata=[1001*201*251]=[T*k1*loglambda])
	! 
	subroutine h5_read3(filenamerd, dataset, rdata)

		IMPLICIT NONE
		CHARACTER(LEN=*)            	:: filenamerd  					! File name for read
		CHARACTER(LEN=*)	     	:: dataset      				! Dataset name
		INTEGER         , PARAMETER 	:: k1dim = 101
		INTEGER         , PARAMETER 	:: k2dim = 501
		INTEGER		, PARAMETER 	:: lambdadim = 136
		INTEGER 			:: error 					! Error flag
		INTEGER(HID_T) 			:: filerd_id, filewr_id 			! file handlers
		INTEGER(HID_T) 			:: dsetrd_id, dsetwr_id, dimspace_id 		! Nint data handles
		INTEGER(HSIZE_T), DIMENSION(1:3):: dims = (/lambdadim, k1dim, k2dim/)	 	! grid data size buffer
		integer 			:: i, j, k					! iteration number
	     	REAL(KIND=8) 	, DIMENSION(lambdadim,k1dim,k2dim), intent(out) :: rdata	! read dataset buffer

		! Initialize FORTRAN interface.
		! 
		CALL h5open_f(error)
		CALL h5fopen_f (filenamerd, H5F_ACC_RDONLY_F, filerd_id, error)

		! Open an existing dataset.
		! 
		CALL h5dopen_f(filerd_id, dataset, dsetrd_id, error)
		CALL h5dread_f(dsetrd_id, H5T_IEEE_F64LE, rdata, dims, error)
		CALL h5dclose_f(dsetrd_id, error)

		! Close the file.
		! 
		CALL h5fclose_f(filerd_id, error)
		!CALL h5close_f(error)

	END SUBROUTINE h5_read3

	!-----------------------------------------------------------------------------------------------------


	! 
	! read all table value for 2D data (rdata2=[1001*201]=[T*k1])
	! 
	subroutine h5_read2(filenamerd, dataset, rdata2)

		IMPLICIT NONE
		CHARACTER(LEN=*), intent(in)   	:: filenamerd  					! File name for read
		CHARACTER(LEN=*), intent(in)  	:: dataset      				! Dataset name
		INTEGER         , PARAMETER 	:: k1dim = 101
		INTEGER         , PARAMETER 	:: k2dim = 501
		INTEGER 			:: error 					! Error flag
		INTEGER(HID_T) 			:: filerd_id, filewr_id 			! file handlers
		INTEGER(HID_T) 			:: dsetrd_id, dsetwr_id, dimspace_id 		! Nint data handles
		INTEGER(HSIZE_T), DIMENSION(1:2):: dims = (/k1dim, k2dim/)	 		! grid data size buffer
		integer 			:: i, j						! iteration number
	     	REAL(KIND=8) 	, DIMENSION(k1dim,k2dim), intent(out) :: rdata2			! read dataset buffer

		! Initialize FORTRAN interface.
		! 
		CALL h5open_f(error)
		CALL h5fopen_f (filenamerd, H5F_ACC_RDONLY_F, filerd_id, error)

		! Open an existing dataset.
		! 
		CALL h5dopen_f(filerd_id, dataset, dsetrd_id, error)
		CALL h5dread_f(dsetrd_id, H5T_IEEE_F64LE, rdata2, dims, error)
		CALL h5dclose_f(dsetrd_id, error)

		! Close the file.
		! 
		CALL h5fclose_f(filerd_id, error)
		!CALL h5close_f(error)

	END SUBROUTINE h5_read2

	!-----------------------------------------------------------------------------------------------------


	!	
	! read all grid data for 3D data (rdata=[1001*201*251]=[T*k1*loglambda]) 
	! 
	subroutine h5_gread3(filenamerd, Tdata, k1data, ldata, Tarray, k1array, loglambdaarray)

		IMPLICIT NONE
		CHARACTER(LEN=*)            	:: filenamerd							! File name for read
		CHARACTER(LEN=*)	     	:: Tdata							! T Dataset name
		CHARACTER(LEN=*)	     	:: k1data							! k1 Dataset name
		CHARACTER(LEN=*)	     	:: ldata 							! loglambda Dataset name
		INTEGER         , PARAMETER 	:: k1dim = 101							! k1 grid size
		INTEGER         , PARAMETER 	:: Tdim = 501							! T grid size
		INTEGER		, PARAMETER 	:: loglambdadim = 136						! loglambda grid size
		INTEGER 			:: error 							! Error flag
		INTEGER(HID_T) 			:: filerd_id					 		! file handles
		INTEGER(HID_T) 			:: dsetT_id = 11, dimspaceT_id = 12			 	! T data handles
		INTEGER(HID_T) 			:: dsetk_id = 21, dimspacek_id = 22			 	! k1 data handles
		INTEGER(HID_T) 			:: dsetl_id = 31, dimspacel_id = 32			 	! loglambda data handles
		INTEGER(HSIZE_T), DIMENSION(1)	:: Tdims = (/Tdim/)	 					! T grid data size buffer
		INTEGER(HSIZE_T), DIMENSION(1)	:: kdims = (/k1dim/)	 					! k1 grid data size buffer
		INTEGER(HSIZE_T), DIMENSION(1)	:: ldims = (/loglambdadim/)	 				! loglambda grid data size buffer
		integer 			:: i, j, k							! iteration number
		REAL(KIND=8) 	, DIMENSION(k1dim), intent(out) :: k1array					! read dataset buffer
		REAL(KIND=8) 	, DIMENSION(Tdim), intent(out) :: Tarray					! read dataset buffer
		REAL(KIND=8) 	, DIMENSION(loglambdadim), intent(out) :: loglambdaarray			! read dataset buffer

		! Initialize FORTRAN interface.
		! 
		CALL h5open_f(error)
		CALL h5fopen_f (filenamerd, H5F_ACC_RDONLY_F, filerd_id, error)

		! Read T array
		! 
		CALL h5dopen_f(filerd_id, Tdata, dsetT_id, error)
		CALL h5dread_f(dsetT_id, H5T_IEEE_F64LE, Tarray, Tdims, error)
		CALL h5dclose_f(dsetT_id, error)

		! Read k1 array
		! 
		CALL h5dopen_f(filerd_id, k1data, dsetk_id, error)
		CALL h5dread_f(dsetk_id, H5T_IEEE_F64LE, k1array, kdims, error)
		CALL h5dclose_f(dsetk_id, error)

		! Read loglambda array
		! 
		CALL h5dopen_f(filerd_id, ldata, dsetl_id, error)
		CALL h5dread_f(dsetl_id, H5T_IEEE_F64LE, loglambdaarray, ldims, error)
		CALL h5dclose_f(dsetl_id, error)
 
		! Close the file. and FORTRAN interface
		! 
		CALL h5fclose_f(filerd_id, error)
		!CALL h5close_f(error)

	END SUBROUTINE h5_gread3

	!-----------------------------------------------------------------------------------------------------

	!	
	! read all grid data for 2D data (Tarray=[1*1001]  k1array=[1*201]) 
	! 
	subroutine h5_gread2(filenamerd, Tdata, k1data, Tarray, k1array)

		IMPLICIT NONE
		CHARACTER(LEN=*)            	:: filenamerd							! File name for read
		CHARACTER(LEN=*)	     	:: Tdata							! T Dataset name
		CHARACTER(LEN=*)	     	:: k1data							! k1 Dataset name
		INTEGER         , PARAMETER 	:: k1dim = 101							! k1 grid size
		INTEGER         , PARAMETER 	:: Tdim = 501							! T grid size
		INTEGER 			:: error 							! Error flag
		INTEGER(HID_T) 			:: filerd_id					 		! file handles
		INTEGER(HID_T) 			:: dsetT_id = 11, dimspaceT_id = 12			 	! T data handles
		INTEGER(HID_T) 			:: dsetk_id = 21, dimspacek_id = 22			 	! k1 data handles
		INTEGER(HSIZE_T), DIMENSION(1)	:: Tdims = (/Tdim/)	 					! T grid data size buffer
		INTEGER(HSIZE_T), DIMENSION(1)	:: kdims = (/k1dim/)	 					! k1 grid data size buffer
		REAL(KIND=8) 	, DIMENSION(k1dim), intent(out) :: k1array					! read dataset buffer
		REAL(KIND=8) 	, DIMENSION(Tdim), intent(out) :: Tarray					! read dataset buffer

		! Initialize FORTRAN interface.
		! 
		CALL h5open_f(error)
		CALL h5fopen_f (filenamerd, H5F_ACC_RDONLY_F, filerd_id, error)
		
		! Read T array
		! 
		CALL h5dopen_f(filerd_id, Tdata, dsetT_id, error)
		CALL h5dread_f(dsetT_id, H5T_IEEE_F64LE, Tarray, Tdims, error)
		CALL h5dclose_f(dsetT_id, error)

		! Read k1 array
		! 
		CALL h5dopen_f(filerd_id, k1data, dsetk_id, error)
		CALL h5dread_f(dsetk_id, H5T_IEEE_F64LE, k1array, kdims, error)
		CALL h5dclose_f(dsetk_id, error)

 
		! Close the file. and FORTRAN interface
		! 
		CALL h5fclose_f(filerd_id, error)
		CALL h5close_f(error)

	END SUBROUTINE h5_gread2


	!-----------------------------------------------------------------------------------------------------	

	! 
	! read T direction spline coefficients for 3D data
	! 
subroutine h5_cread3(cdata, cdataset, bdata, bdataset, ddata, ddataset)

		IMPLICIT NONE
		CHARACTER(LEN=21)            	:: filenamerd = "interpcoefficients.h5" 	! File name for read
		CHARACTER(LEN=*)	     	:: cdataset		  			! Dataset name
		CHARACTER(LEN=*)	     	:: bdataset					! Dataset name
		CHARACTER(LEN=*)	     	:: ddataset					! Dataset name
		INTEGER         , PARAMETER 	:: k1dim = 101
		INTEGER         , PARAMETER 	:: k2dim = 501
		INTEGER		, PARAMETER 	:: lambdadim = 136
		INTEGER 			:: error 					! Error flag
		INTEGER(HID_T) 			:: filerd_id 					! file handlers
		INTEGER(HID_T) 			:: dsetc_id, dimspacec_id		! Nint data handles
		INTEGER(HID_T) 			:: dsetb_id, dimspaceb_id		! Nint data handles
		INTEGER(HID_T) 			:: dsetd_id, dimspaced_id		! Nint data handles
		INTEGER(HSIZE_T), DIMENSION(1:3):: cdims = (/lambdadim, k1dim, k2dim/)	 	! grid data size buffer
		INTEGER(HSIZE_T), DIMENSION(1:3):: bdims = (/lambdadim, k1dim, k2dim-1/)	! grid data size buffer
		INTEGER(HSIZE_T), DIMENSION(1:3):: ddims = (/lambdadim, k1dim, k2dim-1/)	! grid data size buffer
		integer 			:: i, j, k					! iteration number
	     	REAL(KIND=8) 	, DIMENSION(lambdadim,k1dim,k2dim), intent(out) :: cdata	! read dataset buffer
	     	REAL(KIND=8) 	, DIMENSION(lambdadim,k1dim,k2dim-1), intent(out) :: bdata	! read dataset buffer
	     	REAL(KIND=8) 	, DIMENSION(lambdadim,k1dim,k2dim-1), intent(out) :: ddata	! read dataset buffer
		! Initialize FORTRAN interface and open fie
		! 
		CALL h5open_f(error)
		CALL h5fopen_f (filenamerd, H5F_ACC_RDONLY_F, filerd_id, error)
 
		! read c coefficient
		! 
		CALL h5dopen_f(filerd_id, cdataset, dsetc_id, error)
		CALL h5dread_f(dsetc_id, H5T_IEEE_F64LE, cdata, cdims, error)
		CALL h5dclose_f(dsetc_id, error)

		! read b coefficient
		! 
		CALL h5dopen_f(filerd_id, bdataset, dsetb_id, error)
		CALL h5dread_f(dsetb_id, H5T_IEEE_F64LE, bdata, bdims, error)
		CALL h5dclose_f(dsetb_id, error)

		! read d coefficient
		! 
		CALL h5dopen_f(filerd_id, ddataset, dsetd_id, error)
		CALL h5dread_f(dsetd_id, H5T_IEEE_F64LE, ddata, ddims, error)
		CALL h5dclose_f(dsetd_id, error)

		! Close the file and fortran interface
		! 
		CALL h5fclose_f(filerd_id, error)
		!CALL h5close_f(error)

	END SUBROUTINE h5_cread3

	!-----------------------------------------------------------------------------------------------------


	! 
	! read T direction spline coefficients for 2D data
	! 
subroutine h5_cread2(cdata, cdataset, bdata, bdataset, ddata, ddataset)

		IMPLICIT NONE
		CHARACTER(LEN=21)            	:: filenamerd = "interpcoefficients.h5" 	! File name for read
		CHARACTER(LEN=*)	     	:: cdataset		  			! Dataset name
		CHARACTER(LEN=*)	     	:: bdataset					! Dataset name
		CHARACTER(LEN=*)	     	:: ddataset					! Dataset name
		INTEGER         , PARAMETER 	:: k1dim = 101
		INTEGER         , PARAMETER 	:: k2dim = 501
		INTEGER 			:: error 					! Error flag
		INTEGER(HID_T) 			:: filerd_id 					! file handlers
		INTEGER(HID_T) 			:: dsetc_id, dimspacec_id		! Nint data handles
		INTEGER(HID_T) 			:: dsetb_id, dimspaceb_id		! Nint data handles
		INTEGER(HID_T) 			:: dsetd_id, dimspaced_id		! Nint data handles
		INTEGER(HSIZE_T), DIMENSION(1:2):: cdims = (/k1dim, k2dim/)	 		! grid data size buffer
		INTEGER(HSIZE_T), DIMENSION(1:2):: bdims = (/k1dim, k2dim-1/)			! grid data size buffer
		INTEGER(HSIZE_T), DIMENSION(1:2):: ddims = (/k1dim, k2dim-1/)			! grid data size buffer
		integer 			:: i, j, k					! iteration number
	     	REAL(KIND=8) 	, DIMENSION(k1dim,k2dim), intent(out) :: cdata			! read dataset buffer
	     	REAL(KIND=8) 	, DIMENSION(k1dim,k2dim-1), intent(out) :: bdata		! read dataset buffer
	     	REAL(KIND=8) 	, DIMENSION(k1dim,k2dim-1), intent(out) :: ddata		! read dataset buffer
		! Initialize FORTRAN interface and open fie
		! 
		CALL h5open_f(error)
		CALL h5fopen_f (filenamerd, H5F_ACC_RDONLY_F, filerd_id, error)
 
		! read c coefficient
		! 
		CALL h5dopen_f(filerd_id, cdataset, dsetc_id, error)
		CALL h5dread_f(dsetc_id, H5T_IEEE_F64LE, cdata, cdims, error)
		CALL h5dclose_f(dsetc_id, error)

		! read b coefficient
		! 
		CALL h5dopen_f(filerd_id, bdataset, dsetb_id, error)
		CALL h5dread_f(dsetb_id, H5T_IEEE_F64LE, bdata, bdims, error)
		CALL h5dclose_f(dsetb_id, error)

		! read d coefficient
		! 
		CALL h5dopen_f(filerd_id, ddataset, dsetd_id, error)
		CALL h5dread_f(dsetd_id, H5T_IEEE_F64LE, ddata, ddims, error)
		CALL h5dclose_f(dsetd_id, error)

		! Close the file and fortran interface
		! 
		CALL h5fclose_f(filerd_id, error)
		!CALL h5close_f(error)

	END SUBROUTINE h5_cread2

	!-----------------------------------------------------------------------------------------------------
	subroutine h5_checkread(filenamerd, dataset, checkdata)

		IMPLICIT NONE
		CHARACTER(LEN=*)            	:: filenamerd  					! File name for read
		CHARACTER(LEN=*)	     	:: dataset      				! Dataset name
		INTEGER		, PARAMETER	:: checkdim = 140
		INTEGER 			:: error 					! Error flag
		INTEGER(HID_T) 			:: filerd_id, filewr_id 			! file handlers
		INTEGER(HID_T) 			:: dsetrd_id, dsetwr_id, dimspace_id 		! Nint data handles
		INTEGER(HSIZE_T), DIMENSION(1)  :: checkdims = (/checkdim/)			! check data size buffer
		REAL(KIND=8)    , DIMENSION(checkdim)		  , intent(out) :: checkdata	! read check data buffer

		! Initialize FORTRAN interface.
		! 
		CALL h5open_f(error)
		CALL h5fopen_f (filenamerd, H5F_ACC_RDONLY_F, filerd_id, error)

		! Open an existing dataset.
		! 
		CALL h5dopen_f(filerd_id, dataset, dsetrd_id, error)
		CALL h5dread_f(dsetrd_id, H5T_IEEE_F64LE, checkdata, checkdims, error)
		CALL h5dclose_f(dsetrd_id, error)

		! Close the file.
		! 
		CALL h5fclose_f(filerd_id, error)
		CALL h5close_f(error)

	END SUBROUTINE h5_checkread

	!--------------------------------------------------------------------------------------------------------------
	! write output data in HDF5 file 
	subroutine h5_checkwrite(filenamewr, dataset, wdata, ite)

		IMPLICIT NONE
		CHARACTER(LEN=*) 	     	:: filenamewr  	 				! File name for Write
		CHARACTER(LEN=*)	     	:: dataset      				! Dataset name
		INTEGER				:: ite						! ite=1 for first action, ite=0 for second to last action
												! if ite = 1 : create new file , ite = 0 : open exist file
		!INTEGER         , PARAMETER 	:: rltdim = 140					! dimension of result dataset
		INTEGER 			:: error 					! Error flag
		INTEGER(HID_T) 			:: filerd_id, filewr_id 			! file handlers
		INTEGER(HID_T) 			:: dsetrd_id, dsetwr_id, dimspace_id 		! Nint data handles
		INTEGER(HSIZE_T), DIMENSION(3)  :: rltdims	 				! size buffer
		REAL(KIND=8)    , DIMENSION(3), intent(in) :: wdata				! write dataset buffer
		!Real (kind=8), dimension(:,:,:), allocatable, intent(in) :: wdata

		CALL h5open_f(error)
		if (ite .EQ. 1) then			
			Call h5fcreate_f(filenamewr, H5F_ACC_TRUNC_F, filewr_id, error)
		else
			Call h5fopen_f(filenamewr, H5F_ACC_RDWR_F, filewr_id, error)
		end if

		! Create dataspace and write interp result dataset
		! 
		CALL h5screate_simple_f(3, rltdims, dimspace_id, error)
		CALL h5dcreate_f(filewr_id, dataset, H5T_IEEE_F64LE, dimspace_id, dsetwr_id, error)
		CALL h5dwrite_f(dsetwr_id, H5T_IEEE_F64LE, wdata, rltdims, error)
		CALL h5dclose_f(dsetwr_id, error)
		CALL h5sclose_f(dimspace_id, error)

		! Close the file and Fortran interface
		! 
		CALL h5fclose_f(filewr_id, error)
		!CALL h5close_f(error)

	END subroutine h5_checkwrite

	!-----------------------------------------------------------------------------------------------------------
	!write coefficient of 3D grid to interpcoefficients.h5 HDF5 file
	subroutine h5_cwrite3(cdata, cdataset, bdata, bdataset, ddata, ddataset, ite)

		IMPLICIT NONE
		CHARACTER(LEN=21) 	     	:: filenamewr = "interpcoefficients.h5" 	! File name for Write
		CHARACTER(LEN=*)	     	:: cdataset 		   			! Dataset name
		CHARACTER(LEN=*)	     	:: bdataset     	 			! Dataset name
		CHARACTER(LEN=*)	     	:: ddataset     				! Dataset name
		INTEGER				:: ite						! ite=1 for first action, ite=0 for second to last action
												! if ite = 1 : create new file , ite = 0 : open exist file
		INTEGER         , PARAMETER 	:: k1dim = 101
		INTEGER         , PARAMETER 	:: k2dim = 501					! number of T grid
		INTEGER		, PARAMETER 	:: lambdadim = 136
		INTEGER 			:: error 					! Error flag
		INTEGER(HID_T)  		:: filewr_id 			! file handlers
		INTEGER(HID_T) 			:: dsetc_id,  dimspacec_id	 	! Nint data handles
		INTEGER(HID_T) 			:: dsetb_id,  dimspaceb_id		! Nint data handles
		INTEGER(HID_T) 			:: dsetd_id,  dimspaced_id		! Nint data handles
		INTEGER(HSIZE_T), DIMENSION(3):: cdims = (/k2dim, k1dim, lambdadim/)	 	! size buffer
		INTEGER(HSIZE_T), DIMENSION(3):: bdims = (/k2dim-1, k1dim, lambdadim/)	! size buffer
		INTEGER(HSIZE_T), DIMENSION(3):: ddims = (/k2dim-1, k1dim, lambdadim/)	! size buffer
		REAL(KIND=8)    , DIMENSION(k2dim, k1dim, lambdadim), intent(in)   :: cdata	! write dataset buffer
		REAL(KIND=8)    , DIMENSION(k2dim-1, k1dim, lambdadim), intent(in) :: bdata	! write dataset buffer
		REAL(KIND=8)    , DIMENSION(k2dim-1, k1dim, lambdadim), intent(in) :: ddata	! write dataset buffer


		! Initialize FORTRAN interface and create file 
		! 
		CALL h5open_f(error)
		if (ite .EQ. 1) then			
			Call h5fcreate_f(filenamewr, H5F_ACC_TRUNC_F, filewr_id, error)
		else
			Call h5fopen_f(filenamewr, H5F_ACC_RDWR_F, filewr_id, error)
		end if

		CALL h5screate_simple_f(3, cdims, dimspacec_id, error)
		CALL h5screate_simple_f(3, bdims, dimspaceb_id, error)

		
		CALL h5dcreate_f(filewr_id, cdataset, H5T_NATIVE_DOUBLE, dimspacec_id, dsetc_id, error)
		CALL h5dcreate_f(filewr_id, bdataset, H5T_NATIVE_DOUBLE, dimspaceb_id, dsetb_id, error)
		CALL h5dcreate_f(filewr_id, ddataset, H5T_NATIVE_DOUBLE, dimspaceb_id, dsetd_id, error)		

		CALL h5dwrite_f(dsetc_id, H5T_NATIVE_DOUBLE, cdata, cdims, error)
		CALL h5dwrite_f(dsetb_id, H5T_NATIVE_DOUBLE, bdata, bdims, error)
		CALL h5dwrite_f(dsetd_id, H5T_NATIVE_DOUBLE, ddata, bdims, error)

		CALL h5dclose_f(dsetc_id, error)
		CALL h5dclose_f(dsetb_id, error)
		CALL h5dclose_f(dsetd_id, error)
		

		CALL h5sclose_f(dimspacec_id, error)
		CALL h5sclose_f(dimspaceb_id, error)	


		! Close the file and Fortran interface
		! 
		CALL h5fclose_f(filewr_id, error)
		CALL h5close_f(error)

	END subroutine h5_cwrite3

	!-----------------------------------------------------------------------------------------------------------
	! write coefficient of 3D grid to interpcoefficients.h5 HDF5 file
	subroutine h5_cwrite2(cdata, cdataset, bdata, bdataset, ddata, ddataset, ite)

		IMPLICIT NONE
		CHARACTER(LEN=21) 	     	:: filenamewr = "interpcoefficients.h5" 	! File name for Write
		CHARACTER(LEN=*)	     	:: cdataset	 	   			! C Dataset name
		CHARACTER(LEN=*)	     	:: bdataset					! B Dataset name
		CHARACTER(LEN=*)	     	:: ddataset					! D Dataset name
		INTEGER				:: ite
		INTEGER         , PARAMETER 	:: k1dim = 101
		INTEGER         , PARAMETER 	:: k2dim = 501
		INTEGER 			:: error 					! Error flag
		INTEGER(HID_T) 			:: filewr_id 				! file handlers
		INTEGER(HID_T) 			:: dsetc_id,  dimspacec_id	 	! Nint data handles
		INTEGER(HID_T) 			:: dsetb_id,  dimspaceb_id		! Nint data handles
		INTEGER(HID_T) 			:: dsetd_id,  dimspaced_id		! Nint data handles
		INTEGER(HSIZE_T), DIMENSION(2):: cdims = (/k2dim, k1dim/)	 		! size buffer
		INTEGER(HSIZE_T), DIMENSION(2):: bdims = (/k2dim-1, k1dim/)			! size buffer
		INTEGER(HSIZE_T), DIMENSION(2):: ddims = (/k2dim-1, k1dim/)			! size buffer
		REAL(KIND=8)    , DIMENSION(k2dim, k1dim), intent(in)   :: cdata		! write dataset buffer
		REAL(KIND=8)    , DIMENSION(k2dim-1, k1dim), intent(in) :: bdata		! write dataset buffer
		REAL(KIND=8)    , DIMENSION(k2dim-1, k1dim), intent(in) :: ddata		! write dataset buffer

		! Initialize FORTRAN interface and create file 
		! 
		CALL h5open_f(error)
		if (ite .EQ. 1) then			
			Call h5fcreate_f(filenamewr, H5F_ACC_TRUNC_F, filewr_id, error)
		else
			Call h5fopen_f(filenamewr, H5F_ACC_RDWR_F, filewr_id, error)
		end if



		CALL h5screate_simple_f(2, cdims, dimspacec_id, error)
		CALL h5screate_simple_f(2, bdims, dimspaceb_id, error)

		
		CALL h5dcreate_f(filewr_id, cdataset, H5T_NATIVE_DOUBLE, dimspacec_id, dsetc_id, error)
		CALL h5dcreate_f(filewr_id, bdataset, H5T_NATIVE_DOUBLE, dimspaceb_id, dsetb_id, error)
		CALL h5dcreate_f(filewr_id, ddataset, H5T_NATIVE_DOUBLE, dimspaceb_id, dsetd_id, error)		

		CALL h5dwrite_f(dsetc_id, H5T_NATIVE_DOUBLE, cdata, cdims, error)
		CALL h5dwrite_f(dsetb_id, H5T_NATIVE_DOUBLE, bdata, bdims, error)
		CALL h5dwrite_f(dsetd_id, H5T_NATIVE_DOUBLE, ddata, bdims, error)

		CALL h5dclose_f(dsetc_id, error)
		CALL h5dclose_f(dsetb_id, error)
		CALL h5dclose_f(dsetd_id, error)
		

		CALL h5sclose_f(dimspacec_id, error)
		CALL h5sclose_f(dimspaceb_id, error)	


		! Close the file and Fortran innerface
		! 
		CALL h5fclose_f(filewr_id, error)
		CALL h5close_f(error)



  END subroutine h5_cwrite2


	subroutine h5_write_M3(cdata, cdataset, ite)

		IMPLICIT NONE
		CHARACTER(LEN=15) 	     	:: filenamewr = "M_ab_Compare.h5" 			! File name for Write
		CHARACTER(LEN=*)	     	:: cdataset 		   			! Dataset name

		INTEGER				:: ite						! ite=1 for first action, ite=0 for second to last action
												! if ite = 1 : create new file , ite = 0 : open exist file
		INTEGER         , PARAMETER 	:: k1dim = 1280
		!INTEGER         , PARAMETER 	:: k2dim = 1					! number of T grid
		INTEGER		, PARAMETER 	:: lambdadim = 1280
		INTEGER 			:: error 					! Error flag
		INTEGER(HID_T)  		:: filewr_id 					! file handlers
		INTEGER(HID_T) 			:: dsetc_id,  dimspacec_id	 		! Nint data handles
		!INTEGER(HSIZE_T), DIMENSION(3):: cdims = (/k2dim, k1dim, lambdadim/)		! size buffer
		INTEGER(HSIZE_T), DIMENSION(2):: cdims = (/k1dim, lambdadim/)	 		! size buffer
		REAL(KIND=8)    , DIMENSION( k1dim, lambdadim), intent(in)   :: cdata		! write dataset buffer

		! Initialize FORTRAN interface and create file 
		! 
		CALL h5open_f(error)
		if (ite .EQ. 1) then			
			Call h5fcreate_f(filenamewr, H5F_ACC_TRUNC_F, filewr_id, error)
		else
			Call h5fopen_f(filenamewr, H5F_ACC_RDWR_F, filewr_id, error)
		end if

		CALL h5screate_simple_f(2, cdims, dimspacec_id, error)
		CALL h5dcreate_f(filewr_id, cdataset, H5T_NATIVE_DOUBLE, dimspacec_id, dsetc_id, error)
		CALL h5dwrite_f(dsetc_id, H5T_NATIVE_DOUBLE, cdata, cdims, error)
		CALL h5dclose_f(dsetc_id, error)
		CALL h5sclose_f(dimspacec_id, error)
		! Close the file and Fortran interface
		! 
		CALL h5fclose_f(filewr_id, error)
		!CALL h5close_f(error)

	END subroutine h5_write_M3

	subroutine h5_write_F3(cdata, cdataset, Dim1, Dim2, Dim3, ite)

		IMPLICIT NONE
		CHARACTER(LEN=9) 	     	:: filenamewr = "Figure.h5" 			! File name for Write
		CHARACTER(LEN=*)	     	:: cdataset 		   			! Dataset name
		INTEGER				:: Dim1, Dim2, Dim3
		INTEGER				:: ite						! ite=1 for first action, ite=0 for second to last action
												! if ite = 1 : create new file , ite = 0 : open exist file
		!INTEGER         , PARAMETER 	:: k1dim = 1280
		!INTEGER         , PARAMETER 	:: k2dim = 1					! number of T grid
		!INTEGER	 , PARAMETER 	:: lambdadim = 1280
		INTEGER 			:: error 					! Error flag
		INTEGER(HID_T)  		:: filewr_id 					! file handlers
		INTEGER(HID_T) 			:: dsetc_id,  dimspacec_id	 		! Nint data handles
		!INTEGER(HSIZE_T), DIMENSION(3):: cdims = (/k2dim, k1dim, lambdadim/)		! size buffer
		INTEGER(HSIZE_T), DIMENSION(3):: cdims 	 		! size buffer
		REAL(KIND=8)    , DIMENSION(Dim1, Dim2, Dim3), intent(in)   :: cdata		! write dataset buffer
		cdims= (/Dim1, Dim2, Dim3/)
		! Initialize FORTRAN interface and create file 
		! 
		CALL h5open_f(error)
		if (ite .EQ. 1) then			
			Call h5fcreate_f(filenamewr, H5F_ACC_TRUNC_F, filewr_id, error)
		else
			Call h5fopen_f(filenamewr, H5F_ACC_RDWR_F, filewr_id, error)
		end if

		CALL h5screate_simple_f(3, cdims, dimspacec_id, error)
		CALL h5dcreate_f(filewr_id, cdataset, H5T_NATIVE_DOUBLE, dimspacec_id, dsetc_id, error)
		CALL h5dwrite_f(dsetc_id, H5T_NATIVE_DOUBLE, cdata, cdims, error)
		CALL h5dclose_f(dsetc_id, error)
		CALL h5sclose_f(dimspacec_id, error)
		! Close the file and Fortran interface
		! 
		CALL h5fclose_f(filewr_id, error)
		!CALL h5close_f(error)

	END subroutine h5_write_F3

	subroutine h5_write_F2(cdata, cdataset, Dim1, Dim2, ite)

		IMPLICIT NONE
		CHARACTER(LEN=9) 	     	:: filenamewr = "Figure.h5" 			! File name for Write
		CHARACTER(LEN=*)	     	:: cdataset 		   			! Dataset name
		INTEGER				:: Dim1, Dim2
		INTEGER				:: ite						! ite=1 for first action, ite=0 for second to last action
												! if ite = 1 : create new file , ite = 0 : open exist file
		!INTEGER         , PARAMETER 	:: k1dim = 1280
		!INTEGER         , PARAMETER 	:: k2dim = 1					! number of T grid
		!INTEGER	 , PARAMETER 	:: lambdadim = 1280
		INTEGER 			:: error 					! Error flag
		INTEGER(HID_T)  		:: filewr_id 					! file handlers
		INTEGER(HID_T) 			:: dsetc_id,  dimspacec_id	 		! Nint data handles
		!INTEGER(HSIZE_T), DIMENSION(3):: cdims = (/k2dim, k1dim, lambdadim/)		! size buffer
		INTEGER(HSIZE_T), DIMENSION(2):: cdims 	 		! size buffer
		REAL(KIND=8)    , DIMENSION(Dim1, Dim2), intent(in)   :: cdata		! write dataset buffer
		cdims= (/Dim1, Dim2/)
		! Initialize FORTRAN interface and create file 
		! 
		CALL h5open_f(error)
		if (ite .EQ. 1) then			
			Call h5fcreate_f(filenamewr, H5F_ACC_TRUNC_F, filewr_id, error)
		else
			Call h5fopen_f(filenamewr, H5F_ACC_RDWR_F, filewr_id, error)
		end if

		CALL h5screate_simple_f(2, cdims, dimspacec_id, error)
		CALL h5dcreate_f(filewr_id, cdataset, H5T_NATIVE_DOUBLE, dimspacec_id, dsetc_id, error)
		CALL h5dwrite_f(dsetc_id, H5T_NATIVE_DOUBLE, cdata, cdims, error)
		CALL h5dclose_f(dsetc_id, error)
		CALL h5sclose_f(dimspacec_id, error)
		! Close the file and Fortran interface
		! 
		CALL h5fclose_f(filewr_id, error)
		!CALL h5close_f(error)

	END subroutine h5_write_F2

	subroutine h5_write_F1(cdata, cdataset, Dim1, ite)

		IMPLICIT NONE
		CHARACTER(LEN=9) 	     	:: filenamewr = "Figure.h5" 			! File name for Write
		CHARACTER(LEN=*)	     	:: cdataset 		   			! Dataset name
		INTEGER				:: Dim1
		INTEGER				:: ite						! ite=1 for first action, ite=0 for second to last action
												! if ite = 1 : create new file , ite = 0 : open exist file
		!INTEGER         , PARAMETER 	:: k1dim = 1280
		!INTEGER         , PARAMETER 	:: k2dim = 1					! number of T grid
		!INTEGER	 , PARAMETER 	:: lambdadim = 1280
		INTEGER 			:: error 					! Error flag
		INTEGER(HID_T)  		:: filewr_id 					! file handlers
		INTEGER(HID_T) 			:: dsetc_id,  dimspacec_id	 		! Nint data handles
		!INTEGER(HSIZE_T), DIMENSION(3):: cdims = (/k2dim, k1dim, lambdadim/)		! size buffer
		INTEGER(HSIZE_T), DIMENSION(1):: cdims 	 		! size buffer
		REAL(KIND=8)    , DIMENSION(Dim1), intent(in)   :: cdata		! write dataset buffer
		cdims= (/Dim1/)
		! Initialize FORTRAN interface and create file 
		! 
		CALL h5open_f(error)
		if (ite .EQ. 1) then			
			Call h5fcreate_f(filenamewr, H5F_ACC_TRUNC_F, filewr_id, error)
		else
			Call h5fopen_f(filenamewr, H5F_ACC_RDWR_F, filewr_id, error)
		end if

		CALL h5screate_simple_f(1, cdims, dimspacec_id, error)
		CALL h5dcreate_f(filewr_id, cdataset, H5T_NATIVE_DOUBLE, dimspacec_id, dsetc_id, error)
		CALL h5dwrite_f(dsetc_id, H5T_NATIVE_DOUBLE, cdata, cdims, error)
		CALL h5dclose_f(dsetc_id, error)
		CALL h5sclose_f(dimspacec_id, error)
		! Close the file and Fortran interface
		! 
		CALL h5fclose_f(filewr_id, error)
		!CALL h5close_f(error)

	END subroutine h5_write_F1

	subroutine h5_write_I1(cdata, cdataset, Dim1, ite)

		IMPLICIT NONE
		CHARACTER(LEN=9) 	     	:: filenamewr = "Figure.h5" 			! File name for Write
		CHARACTER(LEN=*)	     	:: cdataset 		   			! Dataset name
		INTEGER				:: Dim1
		INTEGER				:: ite						! ite=1 for first action, ite=0 for second to last action
												! if ite = 1 : create new file , ite = 0 : open exist file
		!INTEGER         , PARAMETER 	:: k1dim = 1280
		!INTEGER         , PARAMETER 	:: k2dim = 1					! number of T grid
		!INTEGER	 , PARAMETER 	:: lambdadim = 1280
		INTEGER 			:: error 					! Error flag
		INTEGER(HID_T)  		:: filewr_id 					! file handlers
		INTEGER(HID_T) 			:: dsetc_id,  dimspacec_id	 		! Nint data handles
		!INTEGER(HSIZE_T), DIMENSION(3):: cdims = (/k2dim, k1dim, lambdadim/)		! size buffer
		INTEGER(HSIZE_T), DIMENSION(1):: cdims 	 		! size buffer
		INTEGER   , 	 intent(in)   :: cdata		! write dataset buffer
		cdims= (/Dim1/)
		! Initialize FORTRAN interface and create file 
		! 
		CALL h5open_f(error)
		if (ite .EQ. 1) then			
			Call h5fcreate_f(filenamewr, H5F_ACC_TRUNC_F, filewr_id, error)
		else
			Call h5fopen_f(filenamewr, H5F_ACC_RDWR_F, filewr_id, error)
		end if

		CALL h5screate_simple_f(1, cdims, dimspacec_id, error)
		CALL h5dcreate_f(filewr_id, cdataset, H5T_STD_U32LE, dimspacec_id, dsetc_id, error)
		CALL h5dwrite_f(dsetc_id, H5T_STD_U32LE, cdata, cdims, error)
		CALL h5dclose_f(dsetc_id, error)
		CALL h5sclose_f(dimspacec_id, error)
		! Close the file and Fortran interface
		! 
		CALL h5fclose_f(filewr_id, error)
		!CALL h5close_f(error)

	END subroutine h5_write_I1

	subroutine h5_write_I2(cdata, cdataset, Dim1, Dim2, ite)

		IMPLICIT NONE
		CHARACTER(LEN=9) 	     	:: filenamewr = "Figure.h5" 			! File name for Write
		CHARACTER(LEN=*)	     	:: cdataset 		   			! Dataset name
		INTEGER				:: Dim1, Dim2
		INTEGER				:: ite						! ite=1 for first action, ite=0 for second to last action
												! if ite = 1 : create new file , ite = 0 : open exist file
		!INTEGER         , PARAMETER 	:: k1dim = 1280
		!INTEGER         , PARAMETER 	:: k2dim = 1					! number of T grid
		!INTEGER	 , PARAMETER 	:: lambdadim = 1280
		INTEGER 			:: error 					! Error flag
		INTEGER(HID_T)  		:: filewr_id 					! file handlers
		INTEGER(HID_T) 			:: dsetc_id,  dimspacec_id	 		! Nint data handles
		!INTEGER(HSIZE_T), DIMENSION(3):: cdims = (/k2dim, k1dim, lambdadim/)		! size buffer
		INTEGER(HSIZE_T), DIMENSION(2):: cdims 	 		! size buffer
		INTEGER    , DIMENSION(Dim1, Dim2), intent(in)   :: cdata		 	! write dataset buffer
		cdims= (/Dim1, Dim2/)
		! Initialize FORTRAN interface and create file 
		! 
		CALL h5open_f(error)
		if (ite .EQ. 1) then			
			Call h5fcreate_f(filenamewr, H5F_ACC_TRUNC_F, filewr_id, error)
		else
			Call h5fopen_f(filenamewr, H5F_ACC_RDWR_F, filewr_id, error)
		end if

		CALL h5screate_simple_f(2, cdims, dimspacec_id, error)
		CALL h5dcreate_f(filewr_id, cdataset, H5T_STD_U32LE, dimspacec_id, dsetc_id, error)
		CALL h5dwrite_f(dsetc_id, H5T_STD_U32LE, cdata, cdims, error)
		CALL h5dclose_f(dsetc_id, error)
		CALL h5sclose_f(dimspacec_id, error)
		! Close the file and Fortran interface
		! 
		CALL h5fclose_f(filewr_id, error)
		!CALL h5close_f(error)

	END subroutine h5_write_I2


	subroutine h5_write_M1(cdata, cdataset, ite)

		IMPLICIT NONE
		CHARACTER(LEN=17) 	     	:: filenamewr = "GaussChebyshev.h5" 			! File name for Write
		CHARACTER(LEN=*)	     	:: cdataset 		   			! Dataset name
		INTEGER				:: ite						! ite=1 for first action, ite=0 for second to last action
												! if ite = 1 : create new file , ite = 0 : open exist file
		INTEGER         , PARAMETER 	:: k1dim = 1000
		!INTEGER         , PARAMETER 	:: k2dim = 1					! number of T grid
		!INTEGER		, PARAMETER 	:: lambdadim = 10
		INTEGER 			:: error 					! Error flag
		INTEGER(HID_T)  		:: filewr_id 					! file handlers
		INTEGER(HID_T) 			:: dsetc_id,  dimspacec_id	 		! Nint data handles
		!INTEGER(HSIZE_T), DIMENSION(3):: cdims = (/k2dim, k1dim, lambdadim/)		! size 	buffer
		!INTEGER(HSIZE_T), DIMENSION(2):: cdims = (/k1dim, lambdadim/)	 		! size buffer
		INTEGER(HSIZE_T), DIMENSION(1):: cdims = (/k1dim/)	 			! size buffer		
		REAL(KIND=8)    , DIMENSION( k1dim), intent(in)   :: cdata			! write dataset buffer

		! Initialize FORTRAN interface and create file 
		! 
		CALL h5open_f(error)
		if (ite .EQ. 1) then			
			Call h5fcreate_f(filenamewr, H5F_ACC_TRUNC_F, filewr_id, error)
		else
			Call h5fopen_f(filenamewr, H5F_ACC_RDWR_F, filewr_id, error)
		end if

		CALL h5screate_simple_f(1, cdims, dimspacec_id, error)
		CALL h5dcreate_f(filewr_id, cdataset, H5T_NATIVE_DOUBLE, dimspacec_id, dsetc_id, error)
		CALL h5dwrite_f(dsetc_id, H5T_NATIVE_DOUBLE, cdata, cdims, error)
		CALL h5dclose_f(dsetc_id, error)
		CALL h5sclose_f(dimspacec_id, error)
		! Close the file and Fortran interface
		! 
		CALL h5fclose_f(filewr_id, error)
		!CALL h5close_f(error)

	END subroutine h5_write_M1



END module H5_RDWT



