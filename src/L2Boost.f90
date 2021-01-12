MODULE variables
 implicit none

! Número de ANIMALES, DATOS Y EFECTOS
   integer:: n_animal,n_datos,ntrait,maxround,n_train,n_tun,fold
   real*8:: mtry
   character(len=2)::arg
   character(len=68), allocatable::id(:),id_test(:)

!********************************
!**********PARAMETROS ***********
!********************************
     PARAMETER (ntrait=1)        !Numero de caracteres
!********************************

   integer, allocatable:: ani (:)
   integer, allocatable:: cod(:,:),cen(:) !maxan x número de efectos
   real*8, allocatable:: y(:),y1(:)

! ECUACIONES
   integer:: neq !número de ecuaciones del MME
   integer,allocatable:: snp(:,:)
   real*8, allocatable:: solF(:),sol (:),mse(:),rhs (:)

!GIBBS SAMPLING
   real*8::y_mean,h_c

!VARIABLES DE INFORMACION GENOMICA
	integer ::                       n_snp
    integer, allocatable::           coord(:,:),sel(:),tun(:)
    real*8,  allocatable::           p(:),yp(:),g(:),itcp(:),beta(:),freq(:)
    real*8,  allocatable::           y_hat(:),f(:),USM(:,:),K(:,:),hsmooth(:),dist(:)
	real*8::                         y_m,x_m,x_2

!TRAINING, TUNING, TESTING sets
    integer::n_test,tun_fold
    integer, allocatable::           snp_test(:,:)
    real*8, allocatable::            g_test(:),y_test(:),mse_test(:)
	real*8::                         min_mse,mse_tun,mse_opt,AIC_c_last,AIC_c,trace

!OTROS (contadores,semillas, formatos, etc)
	integer ::                      io,contador,idum,horas,minutos,s,r,t,i,j,jj,q,df,n_df,round,iter_opt
	character (len=8) ::            fecha
	character (len=6) ::            name,learner
	character (len=34)::            fdatos,fpar,ftesting,file_name
	real*8::                        temp,tiempo1,tiempo2,tiempo
	real*8::                        den1,den0,ratio,sc_g,sc_p,sc_e,sc_r,sc_S
	real*8::                        VES,VEf,vg,ve
    real*8::                        x2,x1,min,max,suma,v


end MODULE variables

MODULE inicial
!*****************************************************************************
!*This module sets the initial values, data files, gibb-sampling parameters, *
!* genomic information, number of SNPs, nature of the traits and output files* 
!*****************************************************************************
USE variables
CONTAINS

SUBROUTINE inic

!NAME OF PARAMETER file
 fpar='parBOOST' 
 open (22,file=fpar, form='formatted', status='old')
!Data file name:
read (22,*)   !Jump line
read (22,'(a34)') fdatos 
print *,'DATA FILE=', fdatos
 open (10,file=fdatos, form='formatted', status='old')
!Testing file name:
read (22,*)   !Jump line
read (22,'(a34)') ftesting 
print *,'TESTING FILE=', ftesting
!Model information: number of total effects, covariates, and random effects
read (22,*)   !Jump line
read (22,*) n_snp
print *,'# SNPs=', n_snp
read (22,*) mtry
print *,'mtry=', mtry
read (22,*) v
print *,'smoothing parameter=', v
!INFORMATION ON THE THE TRAIT TO BE ANALYZED
read (22,*)   !Jump line
read (22,*) name !name
!missing values
   read (22,*) df 
   if (df.eq.1) then
      read (22,*) cat_df !identify value
	  print *,cat_df,' will be consider as a missing value in data file'
   else
     print *,'No missing values'
   end if
read (22,*)   !Jump line
read (22,*) maxround !maximun number of rouns allowed

read (22,*)   !Jump line
read (22,*) learner !type of learner used (NPR, OLS)


close (22)
END SUBROUTINE inic
END MODULE inicial


MODULE npr
   use variables
contains
!____________SE IMPLEMENTA LA PARTE NO PARAMETRICA_________________________
!Kernel approach according Nadaraya (1964) and Watson (1964)
!g(x)= y*p / p

subroutine rkhs
integer::uu,n
real*8:: unif
!hsmooth=tamaño de la ventana del estimador kernel

! building the USM norms for each SNPs (Almeida & Vinga, 2002; BMC Bioinformatics 3:6)
    write (*,*) 'SET UP KERNEL DISTANCES USING UNIVERSAL SEQUENCE MAP.'

   	uu=3; !number of unique haplotype segments (0,1,2)
	n=2; !minimum number of dimensions to accommodate uu unique units ceil(log2(uu))
	allocate (USM(n_snp,0:n))

DO s=1,n_snp
	USM(s,0)=unif(x1)
	do t=1,n
		USM(s,t)=0.5d0*( USM(s,t-1)+unif(x1) )
	enddo
ENDDO

end subroutine rkhs
END MODULE npr

MODULE prediction
use variables
contains
subroutine predict_npr
 !PREDICTIONS

 open (37,file=ftesting, form='formatted',status='old')
   file_name = 'ESTIMATES.test_NPR' // arg
 open (38,file=file_name, form='formatted')
   file_name = 'MSE.test_NPR' // arg
 open (39,file=file_name, form='formatted')
   file_name = 'EGBV_NPR.' // arg
 open (36,file=file_name, form='formatted')

write(*,*)  '**Calculating predictions in the testing set**'
PRINT '(a34,i6)','Number of data in testing file =',n_test
io=0;
round=round-1

DEALLOCATE (p,yp,K) !Genomic information
ALLOCATE (p(n_test),yp(n_test),K(n_test,n_animal)) !Genomic information
g_test=0.d0;g=0.d0;mse_test=0.d0;y_hat=y1
   DO s=1,iter_opt !FOR EACH SNP
   	!USE NADARAYA-WATSON ESTIMATOR TO ESTIMATE g(x) FOR INDIVIDUAL r
	do i=1,n_datos
		p(i)=0.d0
		yp(i)=0.d0
		do t=1,n_datos
			K(i,t)=4.d0-2*abs( snp_test(r,sel(s))-snp(t,sel(s)) )
			p(i)= p(i)+K(i,t)
			yp(i)= yp(i)+y_hat(i)*K(i,t)
		enddo
		g(i)=g(i)+v*yp(i)/p(i)
	end do !over t animal
	do r=1,n_test
		p(r)=0.d0
		yp(r)=0.d0
		do t=1,n_datos
			K(r,t)=4.d0-2*abs( snp_test(r,sel(s))-snp(t,sel(s)) )
			p(r)= p(r)+K(r,t)
			yp(r)= yp(r)+y_hat(t)*K(r,t)
		end do !over t animal
		g_test(r)=g_test(r)+v*yp(r)/p(r)
		mse_test(s)=mse_test(s)+(y_test(r)-v*g_test(r))**2
	end do !Over r animal
	write (39,'(i5,f15.3)') s,mse_test(s)/float(n_test)
	do i=1,n_datos
		y_hat(i)=y_hat(i)-g(i)
	end do !over t animal
   ENDDO

   !WRITE GENOMIC ESTIMATES
   do i=1,n_datos
   	write (36,'(a68,f12.3,1x,f12.3)') id(i),y1(i),g(i) 
   enddo
   do i=1,n_test
   	write (38,*) id_test(i),y_test(i),g_test(i) 
   enddo
end subroutine predict_npr

subroutine predict_OLS

integer:: quitar

 !PREDICTIONS
  file_name = 'ESTIMATES.test_OLS' // arg
 open (38,file=file_name, form='formatted')
   file_name = 'MSE.test_OLS' // arg
 open (39,file=file_name, form='formatted')
   file_name = 'EGBV_OLS.' // arg
 open (36,file=file_name, form='formatted')
 
write(*,*)  '**Calculating predictions in the testing set**'
write(*,*)  '** at iteration (optimal iteration)',iter_opt, '**'
write(*,*) 
write (*,'(a34,i6)') 'Number of data in testing file =',n_test

io=0
!Initialize variables
g_test=0.d0;g=0.d0

 DO s=1,iter_opt !FOR EACH SNP
	do j=1,n_datos
		g(j)=g(j) + v*( itcp(s)+beta(s)*snp(j,sel(s)) )
	 enddo
	mse_test=0.d0;quitar=0
	do j=1,n_test
		g_test(j)=g_test(j) + v*( itcp(s)+beta(s)*snp_test(j,sel(s)) )
		if (y_test(j).eq.-9999) then
			quitar=quitar+1
	    		cycle
		endif
		mse_test(s)=mse_test(s)+(y_test(j)-g_test(j))**2
	 enddo
	 write (39,'(i5,f15.3)') s,mse_test(s)/float(n_test-quitar)
 ENDDO

 DO i=1,n_datos
 	write (36,'(a68,f12.3,1x,f12.3)') id(i),y1(i),g(i) 
 ENDDO
 DO j=1,n_test
 	write (38,'(a68,f12.3,1x,f12.3)') id_test(j),y_test(j),g_test(j) 
 ENDDO
end subroutine predict_OLS

end MODULE prediction


MODULE BOOSTING
  use variables
  use npr
  use prediction
contains
subroutine NPRBOOST
  file_name = 'solutions.NPR' // arg
  open (29,file=file_name, form='formatted')

WRITE(*,*) '    USING NON-PARAMETRIC REGRESSION'
WRITE(*,*) '    CALCULATING BANDWITH PARAMETERS'

call rkhs !set USM coordinates and norm for kernels (Almeida & Vinga, 2002; BMC Bioinformatics 3:6)
print *,'start rounds'
mse_opt=999999.d0
AIC_c_last=999999.d0
sel=0;y_hat=0.d0
WRITE(*,*) 'Number of records in the training set:',n_train
WRITE(*,*) 'Number of records in the  tuning set:',n_tun

!Calculate bandwith parameter for each predictor as 20% of the range of distances within predictor (Cornillon et al., 2008. Ann Stat. submitted)
  write (*,*) 'calculating smoothing paramters'
   do s=1,n_snp
      !write (*,*) 'calculating h for snp',s
	  min=99999.d0
	  max=-99999.d0
       hsmooth(s)=1.d0 !0.3d0*(max-min)
   enddo

!STARTS BOOSTING ITERATIONS
WRITE(*,*) '-*-*-*BOOSTING STARTS*-*-*-'
v=0.01d0
DO round=1,maxround
   call cpu_time(tiempo1)
   min_mse=999999.d0
   mse=0.D0
   write (*,'(a6,i5)') 'ROUND ',round
   DO s=1,n_snp !FOR EACH SNP
     if (unif(x1).gt.mtry*0.01d0) cycle
    !  write (*,'(a4,i6)') 'snp ',s
      g=0.d0
      !USE NADARAYA-WATSON ESTIMATOR TO ESTIMATE g(x) FOR INDIVIDUAL r
	  !Kernel approach according Nadaraya (1964) and Watson (1964)
      !g(x)= y*p / p
	do r=1,n_datos
		p(r)=0.d0
	        yp(r)=0.d0
		do t=1,n_datos
			if (tun(i).eq.tun_fold) cycle
			K(r,t)=4.d0-2*abs( snp(r,s)-snp(t,s) )
			p(r)= p(r)+K(r,t)
			yp(r)= yp(r)+y(t)*K(r,t)
		end do !over t animal
		g(r)=yp(r)/p(r)
	end do !Over r animal
	  !CALCULATE MEAN SQUARED ERROR
	do i=1,n_datos
		if (tun(i).eq.tun_fold) cycle
		mse(s)=mse(s)+(y(i)-v*g(i))**2
	enddo
	mse(s)=mse(s)/float(n_train)
	!CHOOSES THE PREDICTOR WITH THE SMALLEST MEAN SQUARED ERROR AND SET f(x)=g(x)
	if (mse(s).lt.min_mse) then
		min_mse=mse(s)
	     	sel(round)=s
		f=v*g
	endif
   ENDDO !Over predictors

   !CALCULATE AIC_c CRITERION (REF) TO TUNE THE NUMBER OF BOOSTING ITERATIONS (Bühlmann, 2006. Ann Stat. 34(2):559-583)
   p=0.d0
   yp=0.d0
   mse_tun=0.d0
   !CALCULATE MEAN SQUARED ERROR in testing set
   do i=1,n_datos
	if (tun(i).ne.tun_fold) cycle
      	do t=1,n_datos
		if (tun(t).eq.tun_fold) cycle
		K(i,t)=4.d0-2*abs( snp(i,s)-snp(t,s) )
		p(i)= p(i)+K(i,t)
		yp(i)= yp(i)+y(t)*K(i,t)
	end do !over t animalif (tun(i).eq.tun_fold) cycle
	mse_tun=mse_tun+(y(i)-v*yp(i)/p(i))**2
   enddo
   mse_tun=mse_tun/float(n_tun)
   !CHECK OPTIMAL ITEARATION (minimum MSE) 
   if (mse_tun.lt.mse_opt) then
      iter_opt=round
      mse_opt=mse_tun
   endif
   write(*,'(a4,i4,a6,f12.3,a9,f12.3)') 'round',round,' MSE =',min_mse, ';MSE_tun=',mse_tun !
   write(*,'(a6,i5,a12,f12.3,a11,f12.3)') 'round ',round,' mse_train =',min_mse,'; mse_tun =',mse_tun
   write (*,*) 'optimal iteration =',iter_opt,'mse=',mse_opt
   write (29,'(i7,2f12.3)') sel(round),min_mse,mse_tun
   y=y-f
call cpu_time(tiempo2)
tiempo=tiempo2-tiempo1
!horas=int(tiempo/float(3600))
!minutos=mod(tiempo,3600.0d0)/float(60)
write (*,*) 'time per boosting iteration=', tiempo

END DO

77 continue

write(*,*)'FINISH with ',round-1,' SNPs selected'
call predict_npr
write(*,*)  'DONE'

stop

END SUBROUTINE NPRBOOST

subroutine OLSBOOST
file_name = 'solutions.OLS' // arg
  open (29,file=file_name, form='formatted')

WRITE(*,*) '    USING ORDINARY LEAST SQUARE REGRESSION'
WRITE(*,*) 'Number of records in the training set:',n_train
WRITE(*,*) 'Number of records in the  tuning set:',n_tun

print *,'start rounds'
mse_opt=999999.d0
AIC_c_last=999999.d0
sel=0;beta=0.d0;itcp=0.d0
vg=0.d0

!STARTS BOOSTING ITERATIONS
WRITE(*,*) '-*-*-*BOOSTING STARTS*-*-*-'

DO round=1,maxround
   call cpu_time(tiempo1)
   min_mse=999999.d0
   mse=0.D0
   write (*,'(a6,i5)') 'ROUND ',round
   DO s=1,n_snp !FOR EACH SNP
      if (unif(x1).gt.mtry*0.01d0) cycle
      y_m=0.d0;x_m=0.d0;x_2=0.d0
	  suma=0.d0
      do i=1,n_datos
	     if (tun(i).eq.tun_fold) cycle
		  y_m=y_m+y(i)
		  x_m=x_m+snp(i,s)
	      suma=suma+y(i)*snp(i,s)/float(n_train)
		  x_2=x_2+snp(i,s)*snp(i,s)/float(n_train)
	  enddo
      b=(suma-(y_m/float(n_train)))*(x_m/float(n_train))/(x_2-(x_m/float(n_train))**2) !(sum-y_m*x_m)/( x_2- (x_m/float(n_datos))**2 )
	  a=y_m/float(n_train)-b*x_m/float(n_train)

     !CALCULATE MEAN SQUARED ERROR
	do i=1,n_datos
	  g(i)=a+b*snp(i,s)
	  if (tun(i).eq.tun_fold) cycle
      mse(s)=mse(s)+(y(i)-g(i))**2 !if (i.gt.n_train) mse(s)=mse(s)+(y(i)-g(i))**2
    enddo
      mse(s)=mse(s)/float(n_train)
      !CHOOSES THE PREDICTOR WITH THE SMALLEST MEAN SQUARED ERROR AND SET f(x)=g(x)
	  if (mse(s).lt.min_mse) then
         min_mse=mse(s)
	     sel(round)=s
	     f=v*g
		 itcp(round)=a
		 beta(round)=b
      endif
	ENDDO !Over SNP s

   mse_tun=0.d0;sce=0.d0;ve=0.d0
	do j=1,n_datos
		if (tun(j).eq.tun_fold) then
			mse_tun=mse_tun+(y(j)-f(j))**2
		else
			sce=sce+(y(j)-f(j))**2
			ve=ve+(y(j)-f(j))/float(n_datos-n_tun)
		endif
	!  g(j)=g(j)+itcp(round)+beta(round)*snp(j,sel(round))
    enddo
	mse_tun=mse_tun/float(n_tun)

    ve=sce/float(n_datos-n_tun)-ve*ve


   !CHECK OPTIMAL ITEARATION (minimum MSE) 
   if (mse_tun.lt.mse_opt) then
      iter_opt=round
	  mse_opt=mse_tun
   endif
   !Add variance captured
   vg=vg+2.d0*freq(sel(round))*(1.d0-freq(sel(round)))*(v*beta(round))*(v*beta(round))

   write(*,'(a6,i5,a12,f12.3,a11,f12.3)') 'round ',round,' mse_train =',min_mse,'; mse_tun =',mse_tun
   write (*,*) 'total genomic variance captured=',vg
   write (*,*) 'estimated genomic heritability=',vg/(vg+ve)
   write (*,*) 'optimal iteration =',iter_opt,'mse=',mse_opt
   write (29,'(i7,5f16.4)') sel(round),itcp(round),beta(round),mse_tun,vg,ve
   y=y-f
call cpu_time(tiempo2)
tiempo=tiempo2-tiempo1
!horas=int(tiempo/float(3600))
!minutos=mod(tiempo,3600.0d0)/float(60)
write (*,*) 'time per boosting iteration=', tiempo

END DO

77 continue

write(*,*)'FINISH with ',round-1,' SNPs selected'
call predict_OLS
write(*,*)  'DONE'

stop

END SUBROUTINE OLSBOOST
END MODULE BOOSTING



program main_boost
!Bayesian Analysis applied to Animal Models
!PROGRAMA PRINCIPAL

  !Modulos usados:
  use variables
  use inicial
  use BOOSTING
  real*8::unif
  !Semillas
!  call random_seed(x1)
!  call random_seed(x2)
   x1 = 0.7283928517d+10
   x2 = 0.7283928517d+10
  idum=567*7345

!LECTURA DE PARAMETROS INICIALES
!-->ZONA INTERACTIVA
  CALL inic
  call getarg(1,arg)
  read (arg,'(I10)') fold
!file_name = 'm' // arg // '.txt'
!************************LECTURA DE FICHEROS************************************
write (*,*) 'Running fold number ',fold
!write (*,*) file_name

!PHENOTYPIC DATA
DO  !Read the number of lines in the data file
 read (10,*,iostat=io)
 IF(io.ne.0) EXIT
 n_datos=n_datos+1
END DO
PRINT '(a30,a13,a1,i6)','Number of data in file ',fdatos,'=',n_datos
io=0;rewind (10)
n_animal=n_datos

ALLOCATE (ani(n_datos),snp(n_animal,n_snp),coord(n_animal+1,2),&
   y1(n_datos),y(n_datos),y_hat(n_datos),f(n_datos),cen(n_datos),sol(n_datos),&
   solF(n_datos),id(n_datos),hsmooth(n_snp),tun(n_datos),freq(n_snp),&
   rhs(n_datos),mse(10*n_snp),sel(maxround),itcp(maxround),beta(maxround),mse_test(maxround))
ALLOCATE (p(n_animal),yp(n_animal),g(n_animal),K(n_animal,n_animal),dist(n_animal)) !Genomic information
!Initialize variables
gebv=0.d0;solF=0.d0;sol=0.d0;temp=0.d0;freq=0.d0
cen=0;
tun=0;n_tun=0
!Fold from the testing set to use as tuning set 
!****************
 tun_fold=1 !****
!****************
  print *,'****READING PHENOTYPE DATA FILE****'

  DO i=1,n_datos 
    !read (10,'(f8.3,8x,i8,30000i2,3000i2)',iostat=io) y(i),id(i),snp(i,1:n_snp) !
    read (10,*,iostat=io) y(i),id(i),snp(i,1:n_snp) !
	IF (io.ne.0) EXIT
	y1(i)=y(i)
	cen(i)=0
	temp=unif(x2)
	if (temp.gt.(fold*0.10d0) .and. temp.le.((fold+1)*0.10d0)) then
	  tun(i)=1
	  n_tun=n_tun+1
	endif
    do kk=1,n_snp
      freq(kk)=freq(kk)+snp(i,kk)/(2.d0*n_datos)
	enddo
  END DO

n_train=n_datos-n_tun

  io=0
  DO i=1,3
    print *,''
	print '(a7,i2,a13,a13)','Line: ',i,' from data file ',fdatos
    print '(a6,a1,f8.2,1x,a11,i2)',name,'=',y(i), 'censored =',cen(i)
	DO jj=1,5
	       print '(a7,i2,a1,i8)','SNP ',jj,'=',snp(i,jj)
	END DO
  END DO
!PAUSE
CLOSE (10)

 open (37,file=ftesting, form='formatted',status='old')

n_test=0;io=0
DO  !Read the number of lines in the testing file
 read (37,*,iostat=io)
 IF(io.ne.0) EXIT
 n_test=n_test+1
END DO
PRINT '(a34,i6)','Number of data in testing file =',n_test
io=0;rewind (37)

ALLOCATE (snp_test(n_test,n_snp),y_test(n_test),g_test(n_test),id_test(n_test))

!Initialize variables
g_test=0.d0;mse_test=0
  DO i=1,n_test 
    !read (37,'(f8.3,8x,i8,30000i2,3000i2)',iostat=io) y_test(i),id_test(i),snp_test(i,1:n_snp) !,valor(i,1:n_cov),cod(i,n_cov+2:n_efectos)
    read (37,*,iostat=io) y_test(i),id_test(i),snp_test(i,1:n_snp) !,valor(i,1:n_cov),cod(i,n_cov+2:n_efectos)
	IF (io.ne.0) EXIT
  END DO
close(37)

!********************COMIENZA BOOSTING***********************
if (learner.eq.'NPR') then
   call NPRBOOST
else if (learner.eq.'OLS') then
   call OLSBOOST
else 
   write (*,*) 'Regression method not set, or wrong code'
endif

end program main_boost



!!!SUBROUTINES USED ALONG THE PROGRAM. DON'T MODIFY BELOW THIS LINE!!!!

          
!___________SUB LUNIF___________________________
!Generates a pseudo-random integer number from a!
! uniform distribution between 0 and ll. (s1 =seed) !
      function lunif (s1,ll)                    !
      implicit doubleprecision (a-h,o-z)        !
      doubleprecision s1,unif                   !
      lunif=int(unif(s1)*ll)+1                  !
      return                                    ! 
      end                                       !
!_______________________________________________!

!______________SUB UNIF_________________________
!Generates a pseudo-random number from a        !
! unif distribution between 0 and 1. (s1 =seed) !
      function unif (s1)                        !
       implicit doubleprecision (a-h,o-z)        !
       doubleprecision s1,unif                   !
      s1 = mod (s1*16807.d0,21477483647.d0)   !
      unif = s1 / 21477483647.d0               !
   !call random_number(unif)
   !lsol=(int((s1*16807.0d0)/2147483647.0d0))   !
   !s1=s1*16807.0d0-lsol*21477483647.0d0        !
   !unif=s1/2147483647.0d0                      !
      return                                    !
      end                                       !
!_______________________________________________!

!______________SUB XNOR______________________________
!Computes the abcise for any function with probability
!    equals to "prob"                                ! 
      FUNCTION XNOR(prob)                            !
      implicit doubleprecision (a-h,o-z)             !
      if (prob.lt..5) then                           !
      p=prob                                         !
      else if (prob.gt..5) then                      !
      p=1.-prob                                      !
      else                                           !
      xnor=0.                                        !
      end if                                         !
      t=dsqrt (log(1./(p*p)))                        !  
      x=t-(2.515517+t*(.802853+t*.010328))/&         !
          (1.+t**(1.432788+t*(.189269+t*.001308)))   !
	if (p.eq.prob) then                            !
	  xnor=-x                                      !
      else                                           !
        xnor=x                                       !
      end if                                         !
      return                                         !
      END                                            !
!____________________________________________________!

!_______________X_NORMAL________________
!Generates a pseudo-random number from a!
!    normal distribution with           !
!    mean= 0 and sd=1. (x1 =seed)       !
      FUNCTION XNORMAL (x1)             !
      implicit doubleprecision (a-h,o-z)!
      xnormal=xnor(unif(x1))            !
      return                            !
      END                               !
!_______________________________________!

