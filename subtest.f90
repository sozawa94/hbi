program test
implicit none
real(8)::p,q,x
integer::i
p=1.d0
q=2.d0
do i=1,10
call calc(x)
write(*,*) x
p=p+x
end do
stop
contains
  subroutine calc(x)
    implicit none
    !real(8),intent(in)::p,q
    real(8),intent(out)::x
    x=p*q
  end subroutine
end program
