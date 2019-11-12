! Copyright (c) 2019, Josh Bevan
! A simple driver example showing usage of fenwick_tree module.
program fen_driver
use fenwick_tree
implicit none
integer, parameter :: sz = 50
real r(sz), r2(6)
integer :: i, rate(sz), fen(sz), psum(sz), chs(8), ind
character(len=6) :: fmt="(50I4)"

call random_number(r)
rate = int(r*9 + 1) !rand ints from 1-9
!manually compute the prefix sum array for comparison
do i=1,sz
	psum(i) = sum(rate(1:i))
end do
fen = build_fen(rate)
write(*,fmt) (i,i=1,sz), rate, fen, psum

call random_number(r2)
chs = [1, int(r2*psum(ubound(psum,1))), 99999]
!perform binary search within fenwick tree in O(log N) time
do i=1,ubound(chs,1)
	print *, chs(i),bitsearch(fen, chs(i))
end do

end program fen_driver