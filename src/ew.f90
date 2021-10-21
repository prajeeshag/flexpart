! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

real function ew(x,p)

  !****************************************************************
  !SAETTIGUNGSDAMPFDRUCK UEBER WASSER IN PA. X IN KELVIN.
  !NACH DER GOFF-GRATCH-FORMEL.
  !****************************************************************

  implicit none

  real :: x, y, a, p, f_qvsat!, c, d

  ew=0.
  if(x.le.0.) stop 'sorry: t not in [k]'
  ! Formula of Goff and Gratch (after Murray, 1966)
  ! if (x.lt.273.15) then
  ! ! Above ice
  !   a = 273.15/x
  !   y = -20.947031*a - 3.56654*log(a) - 2.01889049/a
  !   ew = 5.75185606E10*exp(y)
  ! else
  ! ! Above water
  !   a = 373.15/x 
  !   y = -18.1972839*a + 5.02808*log(a) - 70242.1852*exp(-26.1205253/a) + &
  !     58.0691913*exp(-8.03945282*a)
  !   ew = 7.95357242E10*exp(y)
  ! endif

  ! ! Formula of Magnus (after Murray, 1966)
  ! if (x.lt.273.15) then
  ! ! Above ice
  !   ew = 6.1078*exp(21.8745584*(x-273.15)/(x-7.66))
  ! else 
  ! ! Above water
  !   ew = 6.1078*exp(17.2693882*(x-273.15)/(x-35.86))
  ! endif

  ! Formula of Buck 1981
  ew = f_qvsat(p,x)

  ! ! Original
  ! y=373.15/x ! changed 373.16 to 373.15
  ! a=-7.90298*(y-1.)
  ! a=a+(5.02808*alog(y)) ! removed 0.43429*
  ! c=(1.-(1./y))*11.344
  ! c=-1.+(10.**c)
  ! c=-1.3816*c/(10.**7)
  ! d=(1.-y)*3.49149
  ! d=-1.+(10.**d)
  ! d=8.1328*d/(10.**3)
  ! y=a+c+d
  ! ew=101324.6*(10.**y)       ! Saettigungsdampfdruck in Pa

end function ew
