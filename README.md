bsolver.F90 은 #include 라는 전처리기 명령어를 통해 포함된 것이라 makefile 에는 포함 X
    해당 코드를 그대로 복사하는 역할을 함.

blacs_mod.F90
    사용하진 않으나 혹시 모르니 냅둠.

/FP4D/
─ makefile
    # makefile
- makefile.depend
    # OBJ_FILES & Dependencies
- makefiles/
    # 각 환경에서 compile 옵션 및 링크.
    - makefile.hnuc
    - makefile.kairos
- module_FP4D
    # module list to compile FP4D in hnuc env.
- module_FP4D_kairos
    # module list to compile FP4D in kairos env.
─ eqmodule/
    "준혁이 형이 만들어줌. eqmodule.F90 에 의해 묶임."
    eqmodule/eqmodule.F90
    eqmodule/cnst-srt.c
    eqmodule/cnst-srt.h
    eqmodule/xth-grid.c
    eqmodule/umfpack-util.c
    eqmodule/RZ-grid.c
    eqmodule/dqag.c
    eqmodule/geqdsk.c
    eqmodule/geom-cir.c
    eqmodule/diag-gen.c
- etc
    "사용하진 않음"
    interpol.f
- fp4d.F90 # main program

-	omp_module.F90  # omp module
    my_mpi.F90 	    # mpi module
    blacs_mod.F90   # 사용하진 않으나 혹시 몰라 냅둠
    # 모듈들
    readInput.F90		readHDF5Input.F90	    readNetCDFInput.F90		writeHDF5Output.F90 
    FP4D_globals.F90	FP4D_init.F90 		    FP4D_timer_lib.F90		FP4D_Post.F90	 	  
    FP4D_math.F90 		FP4D_equilibrium.F90 	f_module.F90			RHS_module.F90      
    quasi_f_module.F90 			 		
    # 이하 collision 관련 모듈
    col_f_module.F90	imp_RHS.F90	
    elliptics.F90		h5_rdwt.F90			    spline_interpolation.F90 
    bsolver.F90     # include 라는 전처리기 명령어를 통해 포함된 것이라 makefile 에는 포함 X

>>>>>>> Add README file
