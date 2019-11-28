program main
implicit none
integer::i
real(8)::data(4000,7)
open(11,file='top3.dat')
do i=1,4000
read(11,*) data(i,:)
write(*,*) data(i,5),data(i,7)
end do
stop
end program
