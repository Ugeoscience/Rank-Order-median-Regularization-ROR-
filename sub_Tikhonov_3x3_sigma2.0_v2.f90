!Copyright 2008-2011 SEISCOPE project, All rights reserved.
!Copyright 2013-2021 SEISCOPEII project, All rights reserved.
!
!Redistribution and use in source and binary forms, with or
!without modification, are permitted provided that the following
!conditions are met:
!
!   Redistributions of source code must retain the above copyright
!   notice, this list of conditions and the following disclaimer.
!   Redistributions in binary form must reproduce the above
!   copyright notice, this list of conditions and the following
!   disclaimer in the documentation and/or other materials provided
!   with the distribution.
!   Neither the name of the SEISCOPE project nor the names of
!   its contributors may be used to endorse or promote products
!   derived from this software without specific prior written permission.
!
!Warranty Disclaimer:
!THIS SOFTWARE IS PROVIDED BY THE SEISCOPE PROJECT AND CONTRIBUTORS
!"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
!FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
!SEISCOPE PROJECT OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
!INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
!BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
!STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
!IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!POSSIBILITY OF SUCH DAMAGE.

!
! Modified Tikhonov regularization routine for edge-preserving FWI
! replace this file with "sub_Tikhonov.f90"
! and recompile TOY2DAC
!
! For 1D windows, directional Tikhonov weightings 
! lambda_x and lambda_z need to be set separately in fwi_input file
! For 2D window types, there are no particular directional
! penalties: either set lambda_x or lambda_z greater than zero. 
! If both set above zero it would cause double penalization.
!

subroutine sub_Tikhonov_fcost(pbdir,inv)
  
  IMPLICIT NONE
#include "common.h"
#include "pbdirect.h"
  include "inversion.h"

  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !INVERSION
  TYPE (inversion) :: inv
  
  real :: fcost_tikho
  integer :: i,i1,i2, i1s, i1e, i2s, i2e, ii, ij, n, nn, nw, vnn, vrep, cs, current, crow
  real :: csum, temp
! sigma 2
  real,allocatable,dimension(:,:) :: Cx, Cz, cols, m, t
  integer, parameter :: nwindow = 9
  real, parameter :: cf(nwindow) = [0.027631, 0.066282, 0.123832, 0.180174, 0.204164, &
                                    0.180174, 0.123832, 0.066282, 0.027631]
  ! half window size (total wondow size is (2*nwh+1) including the center )
  real, dimension(nwindow) :: values, values_sorted, indices, indices_sorted
  real, dimension(nwindow+1) :: coeff
  integer, parameter :: nwh = int( (nwindow-1)/2 )
  
  !STORE FCOST DATA
  inv%fcost_data=inv%fcost/inv%scalingfactor
  
  allocate( Cz(pbdir%n1*pbdir%n2 , nwindow ) )
  allocate( cols(pbdir%n1*pbdir%n2 , nwindow ) )
  allocate( m(pbdir%n1*pbdir%n2, 1) )
  allocate( t(pbdir%n1*pbdir%n2, 1) )
  
  !COMPUTE REGUL IN Z DIRECTION
  fcost_tikho=0.
  
  nw = INT((sqrt(2*REAL(nwh)+1)-1)/2)
  do i=1, inv%npar
     Cz(:,:) = 0.
	 cols(:,:) = 0.
     m(:,:) = 0.
     i2s = 1
	 i2e = pbdir%n2
	 do i2=i2s, i2e
		i1s = inv%ibathy(i2)
		i1e = pbdir%n1
		do i1 = i1s,i1e
		    current = INT( i1 + (i2-1) * pbdir%n1 )
			m( current , 1 ) = inv%model(i1,i2,i)
            coeff(:) = 0
            coeff(1:nwindow) = cf(:)
			! getting values from neighboring elements in vertical direction
			! nn is the number of available neighbors; the boundaries are
			! handled automatically
			nn = 0
			values(:) = 1.0e32
			values_sorted(:) = 0.0
			indices(:) = 0 
			indices_sorted(:) = 0 
			crow = 0
			do ij = (i2-nw), (i2+nw)
				do ii = (i1-nw),(i1+nw)
					vnn = ( ii + (i2-1)*pbdir%n1 )
                    if ( ii<=i1e .AND. ij>=i2s .AND. ij<=i2e ) THEN
						if ( ii>=inv%ibathy(ij) ) THEN
						   nn = nn + 1
						   values(nn) = inv%model(ii,ij,i) + &
							  ( abs( -1 + LOG10(10.+REAL(vnn))) *1.0e-8)
				           indices(nn) = ( ii + (i2-1)*pbdir%n1 )
						   if (vnn==current) THEN
							  crow = 1
						   endif
						endif
					endif
				enddo
			enddo

			!sorting
			do vnn = 1, nn
				n = MINLOC(values,1)
				values_sorted(vnn) = values(n)
				indices_sorted(vnn) = indices(n)
				values(n) = 1.0e32
			enddo
			
			vrep = -1
			temp = MOD(nn, 2)
			if (temp == 0.0) THEN
			   do vnn = nwindow+1, 1, -1
				   if ( vnn > INT((nwindow+1)/2) ) THEN
					   coeff(vnn) = cf(vnn-1)
					   vrep = 0
				   endif
			   enddo
			endif
			
			cs = INT( ((nwh+1) - ((nn)/2 ) )  + vrep )
			do vnn = 1, nn
			   cols(current, vnn ) = INT(indices_sorted(vnn))
			   if ( (INT(current)).NE.(INT( indices_sorted(vnn) )  )   ) THEN
					Cz(current, vnn ) =  (coeff(cs+vnn) + 1.0e-8)
			   else
					Cz(current, vnn) = 0               
			   endif
			enddo
			
			if (crow == 0) THEN
			   nn = nn + 1
			   Cz(current, nn) = 0
			   cols(current, nn ) = current
			endif
			
			! control step to ensure each row sums to 0
			csum = SUM(Cz(current, :))
			Cz(current, :) = Cz(current, :) * 1.0/csum
			do vnn = 1, nn
				if ( (INT(current)) == (INT( indices_sorted(vnn) )  )   ) THEN
					Cz(current, vnn) = -1
				endif
			enddo
			
		enddo
	 enddo
	 
	 ! calculate C*m
	 t(:,:) = 0.0
	 do ii = 1, pbdir%n1*pbdir%n2
		do vrep = 1, nwindow
			if (cols(ii,vrep) > 0) THEN
			   t(ii,1) = t(ii,1) + Cz(ii,vrep) * m(INT(cols(ii,vrep)),1)
			endif
		enddo
		fcost_tikho = fcost_tikho + inv%lambda(i) * t(ii,1)**2
	 enddo
	 
  enddo
  
  inv%fcost=inv%fcost+1e0/(pbdir%h**2)*inv%lambda_z*0.5*fcost_tikho  

  !STORE FCOST REG Z
  inv%fcost_reg=1e0/(pbdir%h**2)*inv%lambda_z*0.5*fcost_tikho
  deallocate(Cz)
  deallocate(cols)
  deallocate(m)
  deallocate(t)
  
  !COMPUTE REGUL IN X DIRECTION
  allocate( Cx(pbdir%n1*pbdir%n2, nwindow) )
  allocate( cols(pbdir%n1*pbdir%n2, nwindow) )
  allocate( m(pbdir%n1*pbdir%n2, 1) )
  allocate( t(pbdir%n1*pbdir%n2, 1) )
  fcost_tikho=0.
  
  nw = INT((sqrt(2*REAL(nwh)+1)-1)/2)
  
  do i=1, inv%npar
     Cx(:,:) = 0.
	 cols(:,:) = 0.
     m(:,:) = 0.
     i2s = 1
	 i2e = pbdir%n2
	 do i2=i2s, i2e
		i1s = inv%ibathy(i2)
		i1e = pbdir%n1
		do i1 = i1s,i1e
            coeff(:) = 0
            coeff(1:nwindow) = cf(:)
		    current = INT( i1 + (i2-1) * pbdir%n1 )
			m( current , 1 ) = inv%model(i1,i2,i)
			! getting values from neighboring elements in vertical direction
			! nn is the number of available neighbors; the boundaries are
			! handled automatically
			nn = 0
			values(:) = 1.0e32
			values_sorted(:) = 0.0
			indices(:) = 0 
			indices_sorted(:) = 0 
			crow = 0
			do ij = (i1-nw) , (i1+nw)
				do ii = (i2-nw) , (i2+nw)
					vnn = ( ij + (ii-1)*pbdir%n1 )
					if ( ii>=i2s .AND. ii<=i2e .AND. ij<=i1e ) THEN
						if  ( ij>=inv%ibathy(ii) ) THEN
							nn = nn + 1
							values(nn) = inv%model(ij,ii,1) + &
							   ( abs( -1 + LOG10(10.+REAL(vnn))) *1.0e-8)
				            indices(nn) = ( ij + (ii-1)*pbdir%n1 )
							if (vnn==current) THEN
								crow = 1
							endif
						endif
					endif
				enddo
			enddo

			!sorting
			do vnn = 1, nn
				n = MINLOC(values,1)
				values_sorted(vnn) = values(n)
				indices_sorted(vnn) = indices(n)
				values(n) = 1.0e32
			enddo

			vrep = -1
			temp = MOD(nn, 2)
			if (temp == 0.0) THEN
			   do vnn = nwindow+1, 1, -1
				   if ( vnn > INT((nwindow+1)/2) ) THEN
					   coeff(vnn) = cf(vnn-1)
					   vrep = 0
				   endif
			   enddo
			endif
			
				
			cs = INT( ((nwh+1) - ((nn)/2 ) )  + vrep )
			do vnn = 1, nn
			   cols(current, vnn ) = INT(indices_sorted(vnn))
			   if ( (INT(current)).NE.(INT( indices_sorted(vnn) )  )   ) THEN
					Cx(current, vnn ) =  (coeff(cs+vnn) + 1.0e-8)
			   else
					Cx(current, vnn) = 0               
			   endif
			enddo
			
			if (crow == 0) THEN
			   nn = nn + 1
			   Cx(current, nn) = 0
			   cols(current, nn ) = current
			endif
			
			csum = SUM(Cx(current, :))
			Cx(current, :) = Cx(current, :) * 1.0/csum
			do vnn = 1, nn
				if ( (INT(current)) == (INT( indices_sorted(vnn) )  )   ) THEN
					Cx(current, vnn) = -1
				endif
			enddo
			
		enddo
	 enddo
	 
	 ! calculate C*m
	 t(:,:) = 0.0
	 do ii = 1, pbdir%n1*pbdir%n2
		do vrep = 1, nwindow
			if (cols(ii,vrep) > 0) THEN
			   t(ii,1) = t(ii,1) + Cx(ii,vrep) * m(INT(cols(ii,vrep)),1)
			endif
		enddo
		fcost_tikho = fcost_tikho + inv%lambda(i) * t(ii,1)**2
	 enddo
	 
  enddo

  inv%fcost=inv%fcost+1e0/(pbdir%h**2)*inv%lambda_x*0.5*fcost_tikho  
  
  !STORE FCOST REG X
  inv%fcost_reg=inv%fcost_reg+&
       1e0/(pbdir%h**2)*inv%lambda_x*0.5*fcost_tikho  

  !NORMALIZE FCOST REG
  inv%fcost_reg=inv%fcost_reg/inv%scalingfactor
  
  
  deallocate(Cx)
  deallocate(cols)
  deallocate(m)
  deallocate(t)

end subroutine sub_Tikhonov_fcost


subroutine sub_Tikhonov_fgrad(pbdir,inv)
  
  IMPLICIT NONE
#include "common.h"
#include "pbdirect.h"
  include "inversion.h"

  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !INVERSION
  TYPE (inversion) :: inv

  real,allocatable,dimension(:,:,:) :: fgrad_tikho_x,fgrad_tikho_z
  real,allocatable,dimension(:,:) :: Cx, Cz, cols, m, t, t2
  integer :: i,i1,i2, i1s, i1e, i2s, i2e, ii, ij, n, nn, nw, vnn, vrep, cs, current, crow
  real :: csum, temp

! sigma 2
  integer, parameter :: nwindow = 9
  real, parameter :: cf(nwindow) = [0.027631, 0.066282, 0.123832, 0.180174, 0.204164, &
                                    0.180174, 0.123832, 0.066282, 0.027631]
  ! half window size (total wondow size is (2*nwh+1) including the center )
  real, dimension(nwindow) :: values, values_sorted, indices, indices_sorted
  real, dimension(nwindow+1) :: coeff 
  integer, parameter :: nwh = int( (nwindow-1)/2 )
  
  !Vertical term
  allocate( fgrad_tikho_z(pbdir%n1,pbdir%n2,inv%npar) )
  allocate( Cz(pbdir%n1*pbdir%n2 , nwindow ) )
  allocate( cols(pbdir%n1*pbdir%n2 , nwindow ) )
  allocate( m(pbdir%n1*pbdir%n2, 1) )
  allocate( t(pbdir%n1*pbdir%n2, 1) )
  allocate( t2(pbdir%n1*pbdir%n2, 1) )
  fgrad_tikho_z(:,:,:)= 0.
  
  nw = INT((sqrt(2*REAL(nwh)+1)-1)/2)
  
  do i=1, inv%npar
     Cz(:,:) = 0.
	 cols(:,:) = 0.
     m(:,:) = 0.
     i2s = 1
	 i2e = pbdir%n2
	 do i2=i2s, i2e
		i1s = inv%ibathy(i2)
		i1e = pbdir%n1
		do i1 = i1s,i1e
            coeff(:) = 0
            coeff(1:nwindow) = cf(:)
		    current = INT( i1 + (i2-1) * pbdir%n1 )
            m( current , 1 ) = inv%model(i1,i2,i)
			! getting values from neighboring elements in vertical direction
			! nn is the number of available neighbors; the boundaries are
			! handled automatically
			nn = 0
			values(:) = 1.0e32
			values_sorted(:) = 0.0
			indices(:) = 0 
			indices_sorted(:) = 0 
			crow = 0
			do ij = (i2-nw), (i2+nw)
				do ii = (i1-nw),(i1+nw)
					vnn = ( ii + (ij-1)*pbdir%n1 )
                    if ( ii<=i1e .AND. ij>=i2s .AND. ij<=i2e ) THEN
						if ( ii>=inv%ibathy(ij) ) THEN
						   nn = nn + 1
						   values(nn) = inv%model(ii,ij,i) + &
							   ( abs( -1 + LOG10(10.+REAL(vnn))) *1.0e-8)
				           indices(nn) = ( ii + (ij-1)*pbdir%n1 )
						   if (vnn==current) THEN
							   crow = 1
						   endif
						endif
					endif
				enddo
			enddo

			!sorting
			do vnn = 1, nn
				n = MINLOC(values,1)
				values_sorted(vnn) = values(n)
				indices_sorted(vnn) = indices(n)
				values(n) = 1.0e32
			enddo

			vrep = -1
			temp = MOD(nn, 2)
			if (temp == 0.0) THEN
			   do vnn = nwindow+1, 1, -1
				   if ( vnn > INT((nwindow+1)/2) ) THEN
					   coeff(vnn) = cf(vnn-1)
					   vrep = 0
				   endif
			   enddo
			endif
			
			cs = INT( ((nwh+1) - ((nn)/2 ) )  + vrep )
			do vnn = 1, nn
			   cols(current, vnn ) = INT(indices_sorted(vnn))
			   if ( (INT(current)).NE.(INT( indices_sorted(vnn) )  )   ) THEN
					Cz(current, vnn ) =  (coeff(cs+vnn) + 1.0e-8)
			   else
					Cz(current, vnn) = 0               
			   endif
			enddo
			
			if (crow == 0) THEN
			   nn = nn + 1
			   Cz(current, nn) = 0
			   cols(current, nn ) = current
			endif
			
			csum = SUM(Cz(current, :))
			Cz(current, :) = Cz(current, :) * 1.0/csum
			do vnn = 1, nn
				if ( (INT(current)) == (INT( indices_sorted(vnn) )  )   ) THEN
					Cz(current, vnn) = -1
				endif
			enddo
			
		enddo
	 enddo
	 
	 !fgrad_tikho_z(:,:,i) = MATMUL( MATMUL( TRANSPOSE(C), C) , m )
	 ! calc Cm
	 t(:,:) = 0.0
	 do ii = 1, pbdir%n1*pbdir%n2
		do vrep = 1, nwindow
			if (cols(ii,vrep) > 0) THEN
			   t(ii,1) = t(ii,1) + Cz(ii,vrep) * m(INT(cols(ii,vrep)),1)
			endif
		enddo
	 enddo
	 
	 ! calc C'(Cm)
	 t2(:,:) = 0.0
	 do ii = 1, pbdir%n1*pbdir%n2
	 	 do vrep = 1, nwindow
			if (cols(ii,vrep) > 0) THEN
			   n = INT(cols(ii,vrep))
			   t2(n,1) = t2(n,1) + Cz(ii,vrep) * t(ii,1)
			endif
		 enddo
	 enddo
	 ! map back to fgrad_tikho
	 do i2=1, pbdir%n2
		i1s = inv%ibathy(i2)
		i1e = pbdir%n1
		do i1 = i1s,i1e
		    current = INT( i1 + (i2-1) * pbdir%n1 )
			fgrad_tikho_z(i1,i2,i) = t2(current,1)
		enddo
	 enddo
	 
  enddo
  
  
  do i=1,inv%npar
     inv%gradient(:,:,i)=inv%gradient(:,:,i)+inv%lambda(i)/(pbdir%h**2)*inv%lambda_z*fgrad_tikho_z(:,:,i)       
  enddo
  
  deallocate(fgrad_tikho_z)
  deallocate(Cz)
  deallocate(cols)
  deallocate(m)
  deallocate(t)
  deallocate(t2)
  
  !Horizontal term
  allocate(fgrad_tikho_x(pbdir%n1,pbdir%n2,inv%npar))
  allocate( Cx(pbdir%n1*pbdir%n2, nwindow) )
  allocate( cols(pbdir%n1*pbdir%n2, nwindow) )
  allocate( m(pbdir%n1*pbdir%n2, 1) )
  allocate( t(pbdir%n1*pbdir%n2, 1) )
  allocate( t2(pbdir%n1*pbdir%n2, 1) )
  fgrad_tikho_x(:,:,:)=0.
  
  nw = INT((sqrt(2*REAL(nwh)+1)-1)/2)
  do i=1, inv%npar
     Cx(:,:) = 0.
	 cols(:,:) = 0.
     m(:,:) = 0.
     i2s = 1
	 i2e = pbdir%n2
	 do i2=i2s, i2e
		i1s = inv%ibathy(i2)
		i1e = pbdir%n1
		do i1 = i1s,i1e
            coeff(:) = 0
            coeff(1:nwindow) = cf(:)
		    current = INT( i1 + (i2-1) * pbdir%n1 )
            m( current , 1 ) = inv%model(i1,i2,i)
			! getting values from neighboring elements in vertical direction
			! nn is the number of available neighbors; the boundaries are
			! handled automatically
			nn = 0
			values(:) = 1.0e32
			values_sorted(:) = 0.0
			indices(:) = 0 
			indices_sorted(:) = 0 
			crow = 0
			do ij = (i1-nw) , (i1+nw)
				do ii = (i2-nw) , (i2+nw)
					vnn = ( ij + (ii-1)*pbdir%n1 )
					if ( ii>=i2s .AND. ii<=i2e .AND. ij<=i1e ) THEN
						if  ( ij>=inv%ibathy(ii) ) THEN
							nn = nn + 1
							values(nn) = inv%model(ij,ii,1) + &
							    ( abs( -1 + LOG10(10+REAL(vnn))) *1.0e-8)
	                        indices(nn) = ( ij + (ii-1)*pbdir%n1 )
							if (vnn==current) THEN
								crow = 1
							endif
						endif
					endif
				enddo
			enddo
			
			!sorting
			do vnn = 1, nn
				n = MINLOC(values,1)
				values_sorted(vnn) = values(n)
				indices_sorted(vnn) = indices(n)
				values(n) = 1.0e32
			enddo
			
			vrep = -1
			temp = MOD(nn, 2)
			if (temp == 0.0) THEN
			   do vnn = nwindow+1, 1, -1
				   if ( vnn > INT((nwindow+1)/2) ) THEN
					   coeff(vnn) = cf(vnn-1)
					   vrep = 0
				   endif
			   enddo
			endif
			
			cs = INT( ((nwh+1) - ((nn)/2 ) )  + vrep )
			do vnn = 1, nn
			   cols(current, vnn ) = INT(indices_sorted(vnn))
			   if ( (INT(current)).NE.(INT( indices_sorted(vnn) )  )   ) THEN
					Cx(current, vnn ) =  (coeff(cs+vnn) + 1.0e-8)
			   else
					Cx(current, vnn) = 0               
			   endif
			enddo
			
			if (crow == 0) THEN
			   nn = nn + 1
			   Cx(current, nn) = 0
			   cols(current, nn ) = current
			endif
			
			csum = SUM(Cx(current, :))
			Cx(current, :) = Cx(current, :) * 1.0/csum
			do vnn = 1, nn
				if ( (INT(current)) == (INT( indices_sorted(vnn) )  )   ) THEN
					Cx(current, vnn) = -1
				endif
			enddo
			
		enddo
	 enddo
	 
	 !fgrad_tikho_x(:,:,i) = MATMUL( MATMUL( TRANSPOSE(C), C) , m )
	 
	 ! calc Cm
	 t(:,:) = 0.0
	 do ii = 1, pbdir%n1*pbdir%n2
		do vrep = 1, nwindow
			if (cols(ii,vrep) > 0) THEN
			   t(ii,1) = t(ii,1) + Cx(ii,vrep) * m(INT(cols(ii,vrep)),1)
			endif
		enddo
	 enddo
	 
	 ! calc C'(Cm)
	 t2(:,:) = 0.0
	 do ii = 1, pbdir%n1*pbdir%n2
	 	 do vrep = 1, nwindow
			if (cols(ii,vrep) > 0) THEN
			   n = INT(cols(ii,vrep))
			   t2(n,1) = t2(n,1) + Cx(ii,vrep) * t(ii,1)
			endif
		 enddo
	 enddo
	 ! map back to fgrad_tikho
	 do i2=1, pbdir%n2
		i1s = inv%ibathy(i2)
		i1e = pbdir%n1
		do i1 = i1s,i1e
		    current = INT( i1 + (i2-1) * pbdir%n1 )
			fgrad_tikho_x(i1,i2,i) = t2(current,1)
		enddo
	 enddo
	 
  enddo
  
  
  do i=1,inv%npar
     inv%gradient(:,:,i)=inv%gradient(:,:,i)+inv%lambda(i)/(pbdir%h**2)*inv%lambda_x*fgrad_tikho_x(:,:,i)
  enddo
  
  deallocate(fgrad_tikho_x)
  deallocate(Cx)
  deallocate(cols)
  deallocate(m)
  deallocate(t)
  deallocate(t2)
  
end subroutine sub_Tikhonov_fgrad

subroutine sub_Tikhonov_Hv(pbdir,inv,Hv,v)
  
  IMPLICIT NONE
#include "common.h"
#include "pbdirect.h"
  include "inversion.h"

  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !INVERSION
  TYPE (inversion) :: inv
  real,dimension(pbdir%n1*pbdir%n2*inv%npar) :: v,Hv
  real,allocatable,dimension(:) :: vect_x,vect_z
  real,allocatable,dimension(:,:) :: Cz, Cx, cols, m, t, t2
  integer :: i,i1,i2, i1s, i1e, i2s, i2e, ii, ij, n, nm, nm2, nn, nw, &
                     vnn, vrep, last, current, currentn, neighbor, neighborn, cs, crow
  ! coefficients sum up to 2

! sigma 2
  integer, parameter :: nwindow = 9
  real, parameter :: cf(nwindow) = [0.027631, 0.066282, 0.123832, 0.180174, 0.204164, &
                                    0.180174, 0.123832, 0.066282, 0.027631]
  ! half window size (total wondow size is (2*nwh+1) including the center )
  real, dimension(nwindow) :: values, values_sorted, indices, indices_sorted
  real, dimension(nwindow+1) :: coeff
  integer, parameter :: nwh = int( (nwindow-1)/2 )
  real :: csum, temp

   !Vertical term
  allocate(vect_z(pbdir%n1*pbdir%n2*inv%npar))
  allocate( Cz(pbdir%n1*pbdir%n2, nwindow) )
  allocate( cols(pbdir%n1*pbdir%n2, nwindow) )
  allocate( m(pbdir%n1*pbdir%n2, 1) )
  allocate( t(pbdir%n1*pbdir%n2, 1) )
  allocate( t2(pbdir%n1*pbdir%n2, 1) )
  vect_z(:) = 0.

  ! total number of parameters    
  nm = pbdir%n1*pbdir%n2*inv%npar
  nm2 = pbdir%n1*pbdir%n2
  nw = INT((sqrt(2*REAL(nwh)+1)-1)/2)
  
  do i=1, inv%npar 
     i2s=1
     i2e=pbdir%n2
	 Cz(:,:) = 0.0
	 cols(:,:) = 0.0
	 m(:,:) = 0.0 
     do i2 = i2s,i2e  
		 i1s = inv%ibathy(i2)
		 i1e = pbdir%n1
		 last = nm
         do i1 = i1s , i1e
			coeff(:) = 0
			coeff(1:nwindow) = cf(:)
			nn = 0
			values(:) = 1.0e32
			values_sorted(:) = 0.0
			indices(:) = 0 
			indices_sorted(:) = 0 
			current = INT( i1+(i2-1)*pbdir%n1 + &
										   (i-1)*pbdir%n1*pbdir%n2)
			currentn = INT( i1 + (i2-1)*pbdir%n1 )
			m( currentn , 1 ) = v(current)
			crow = 0   
			 do ij = (-nw) , (+nw)
                 do ii = (-nw) , (+nw)
                    neighbor = INT( (i1+ii)+(i2+ij-1)*pbdir%n1 + &
					                         (i-1)*pbdir%n1*pbdir%n2 )
                    neighborn = INT( (i1+ii)+(i2+ij-1)*pbdir%n1  )
                    if ( (i1+ii)<=i1e .AND. (i2+ij)>=i2s .AND. &  
					     (i2+ij)<=i2e .AND. neighborn<=nm2 ) THEN
						if ( (i1+ii)>=inv%ibathy(i2+ij) ) THEN
							nn = nn + 1
							values(nn) = v(neighbor) + &
							   ( abs( -1 + LOG10(10+REAL(neighborn))) *1.0e-8)
							indices(nn) = neighborn
							if (currentn == neighborn) THEN
								crow = 1
							endif
						endif
                    endif
                 enddo
             enddo
             
			!sorting
			do vnn = 1,nn
			 n = MINLOC(values,1)
			 values_sorted(vnn) = values(n)
			 indices_sorted(vnn) = indices(n)
			 values(n) = 1.0e32
			enddo
 
			vrep = -1
			temp = MOD(nn, 2)
			if (temp == 0.0) THEN
			   do vnn = nwindow+1, 1, -1
				   if ( vnn > INT((nwindow+1)/2) ) THEN
					   coeff(vnn) = cf(vnn-1)
					   vrep = 0
				   endif
			   enddo
			endif
			

			cs = INT( ((nwh+1) - ((nn)/2 ) )  + vrep )
			do vnn = 1, nn
			   cols(current, vnn ) = INT(indices_sorted(vnn))
			   if ( (INT(current)).NE.(INT( indices_sorted(vnn) )  )   ) THEN
					Cz(current, vnn ) =  (coeff(cs+vnn) + 1.0e-8)
			   else
					Cz(current, vnn) = 0               
			   endif
			enddo
			
			if (crow == 0) THEN
			   nn = nn + 1
			   Cz(current, nn) = 0
			   cols(current, nn ) = current
			endif
			
			csum = SUM(Cz(current, :))
			Cz(current, :) = Cz(current, :) * 1.0/csum
			do vnn = 1, nn
				if ( (INT(current)) == (INT( indices_sorted(vnn) )  )   ) THEN
					Cz(current, vnn) = -1
				endif
			enddo
			 
			 
         enddo
     enddo
     
	 !m = MATMUL( MATMUL( TRANSPOSE(C), C) , m )
	 
	 ! calc Cm
	 t(:,:) = 0.0
	 do ii = 1, pbdir%n1*pbdir%n2
		do vrep = 1, nwindow
			if (cols(ii,vrep) > 0) THEN
			   t(ii,1) = t(ii,1) + Cz(ii,vrep) * m(INT(cols(ii,vrep)),1)
			endif
		enddo
	 enddo
	 
	 ! calc C'(Cm)
	 t2(:,:) = 0.0
	 do ii = 1, pbdir%n1*pbdir%n2
	 	 do vrep = 1, nwindow
			if (cols(ii,vrep) > 0) THEN
			   n = INT(cols(ii,vrep))
			   t2(n,1) = t2(n,1) + Cz(ii,vrep) * t(ii,1)
			endif
		 enddo
	 enddo
	 
	 ! map back to vect
     i2s=1
	 i2e=pbdir%n2
     do i2 = i2s,i2e  
         i1s = 1
         i1e = pbdir%n1
         do i1 = i1s , i1e
             current = INT( i1+(i2-1)*pbdir%n1 + &
                                               (i-1)*pbdir%n1*pbdir%n2)
             currentn = INT( i1+(i2-1)*pbdir%n1 )
			 vect_z(current) = t2(currentn,1)
         enddo
     enddo
	 
  enddo
    
  Hv(:)=Hv(:)+1e0/(pbdir%h**2)*inv%lambda_z*vect_z(:)
  
  deallocate(vect_z)
  deallocate(Cz)
  deallocate(cols)
  deallocate(m)
  deallocate(t)
  deallocate(t2)
  
  !Horizontal term  
  ! total number of parameters    
  nw = INT((sqrt(2*REAL(nwh)+1)-1)/2)
  nm = pbdir%n1*pbdir%n2*inv%npar
  nm2 = pbdir%n1*pbdir%n2
  
  allocate(vect_x(pbdir%n1*pbdir%n2*inv%npar))
  allocate( Cx(pbdir%n1*pbdir%n2, nwindow) )
  allocate( cols(pbdir%n1*pbdir%n2, nwindow) )
  allocate( m(pbdir%n1*pbdir%n2, 1) )
  allocate( t(pbdir%n1*pbdir%n2, 1) )
  allocate( t2(pbdir%n1*pbdir%n2, 1) )
  vect_x(:)=0.
  
  do i=1,inv%npar  
     i2s=1
     i2e=pbdir%n2
	 Cx(:,:) = 0.0
	 cols(:,:) = 0.0
	 m(:,:) = 0.0 
	 do i2 = i2s,i2e
		i1s = inv%ibathy(i2)
		i1e = pbdir%n1  
		last = nm
		do i1 = i1s , i1e
            coeff(:) = 0
            coeff(1:nwindow) = cf(:)
            nn = 0
            values(:) = 1.0e32
            values_sorted(:) = 0.0
            indices(:) = 0 
            indices_sorted(:) = 0 
			current = INT( i1+(i2-1)*pbdir%n1 + (i-1)*pbdir%n1*pbdir%n2)
			currentn = INT( i1+(i2-1)*pbdir%n1)
            m( currentn , 1 ) = v(current)
			crow = 0	
			do ij = (-nw) , nw
				do ii = (-nw),nw
					neighbor = INT( (i1+ij) + (i2-1 + ii)*pbdir%n1 + &
												(i-1)*pbdir%n1*pbdir%n2 )
					neighborn = INT( (i1+ij) + (i2-1 + ii)*pbdir%n1  )
                    if ( (i1+ij)<=i1e .AND. (i2+ii)>=i2s .AND. &
					     (i2+ii)<=i2e .AND. neighborn<=nm2) THEN 
						if ( (i1+ij)>=inv%ibathy(i2+ii) ) THEN
							nn = nn + 1
							values(nn) = v(neighbor) + &
							    ( abs( -1 + LOG10(10+REAL(neighborn))) *1.0e-8)
						    indices(nn) = neighborn
							if (currentn == neighborn) THEN
								crow = 1
							endif
						endif
					endif
				enddo
			enddo
			
		   !sorting
			do vnn = 1,nn
				n = MINLOC(values,1)
				values_sorted(vnn) = values(n)
				indices_sorted(vnn) = indices(n)
				values(n) = 1.0e32
			enddo

			vrep = -1
			temp = MOD(nn, 2)
			if (temp == 0.0) THEN
			   do vnn = nwindow+1, 1, -1
				   if ( vnn > INT((nwindow+1)/2) ) THEN
					   coeff(vnn) = cf(vnn-1)
					   vrep = 0
				   endif
			   enddo
			endif
			
			cs = INT( ((nwh+1) - ((nn)/2 ) )  + vrep )
			do vnn = 1, nn
			   cols(current, vnn ) = INT(indices_sorted(vnn))
			   if ( (INT(current)).NE.(INT( indices_sorted(vnn) )  )   ) THEN
					Cx(current, vnn ) =  (coeff(cs+vnn) + 1.0e-8)
			   else
					Cx(current, vnn) = 0               
			   endif
			enddo
			
			if (crow == 0) THEN
			   nn = nn + 1
			   Cx(current, nn) = 0
			   cols(current, nn ) = current
			endif
			
			csum = SUM(Cx(current, :))
			Cx(current, :) = Cx(current, :) * 1.0/csum
			do vnn = 1, nn
				if ( (INT(current)) == (INT( indices_sorted(vnn) )  )   ) THEN
					Cx(current, vnn) = -1
				endif
			enddo

		enddo
	 enddo
	 
	 ! MATMUL( MATMUL( TRANSPOSE(C), C) , m )
	 ! calc Cm
	 t(:,:) = 0.0
	 do ii = 1, pbdir%n1*pbdir%n2
		do vrep = 1, nwindow
			if (cols(ii,vrep) > 0) THEN
			   t(ii,1) = t(ii,1) + Cx(ii,vrep) * m(INT(cols(ii,vrep)),1)
			endif
		enddo
	 enddo
	 
	 ! calc C'(Cm)
	 t2(:,:) = 0.0
	 do ii = 1, pbdir%n1*pbdir%n2
	 	 do vrep = 1, nwindow
			if (cols(ii,vrep) > 0) THEN
			   n = INT(cols(ii,vrep))
			   t2(n,1) = t2(n,1) + Cx(ii,vrep) * t(ii,1)
			endif
		 enddo
	 enddo
	 
	 ! map back to vect
	 
	 do i2 = i2s,i2e
		i1s = inv%ibathy(i2)
		i1e = pbdir%n1  
		do i1 = i1s , i1e
			current = INT( i1+(i2-1)*pbdir%n1 + (i-1)*pbdir%n1*pbdir%n2)
			currentn = INT( i1+(i2-1)*pbdir%n1)
			vect_x(current) = t2( currentn , 1 )
		enddo
	 enddo

  enddo
  
  Hv(:)=Hv(:)+1e0/(pbdir%h**2)*inv%lambda_x*vect_x(:)  
  
  deallocate(vect_x)
  deallocate(Cx)
  deallocate(cols)
  deallocate(m)
  deallocate(t)
  deallocate(t2)
     
end subroutine sub_Tikhonov_Hv

