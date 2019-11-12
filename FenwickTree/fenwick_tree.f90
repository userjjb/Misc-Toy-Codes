! Copyright (c) 2019, Josh Bevan
! All rights reserved.
! Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
! 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
! 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
! 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
module fenwick_tree
implicit none

contains
! Build a Fenwick Tree aka Binary Indexed Tree from a parent array.
! If desired this could be modified for commutative operations besides
! prefix sums, and for other data types.
pure function build_fen(arr) result(fen)
	integer, dimension(:), intent(in) :: arr
	integer, dimension(size(arr)) :: fen
	integer :: i,j,n

	fen = 0
	n = size(arr)
	do i=1,n
		j = i
		do while (j<=n)
			fen(j) = fen(j) + arr(i)
			j = j + iand(j,-j)
		end do
	end do
end function build_fen

! Perform a binary search over the Fenwick tree, using the implicit tree
! structure instead of the usual bisection. Cost is therefore O(log N).
! NOTE: This only works for monotonic Fenwick trees.
pure function bitsearch(tree, tgt) result(best)
	integer, intent(in) :: tgt
	integer, dimension(:), intent(in) :: tree
	integer :: ind, best, bit, accum, n

	n = size(tree)
	bit = int(log(n*1.)/log(2.)) !total bits needed-2
	ind = 2**bit !Index = top bit of tree set
	accum = 0; best = 1;
	do while (bit>=0)
		bit = bit-1
		accum = accum + tree(ind)
		if (accum>tgt) then
			accum = accum - tree(ind)
			ind = ind - 2**bit
		else if (accum<tgt) then
			best = ind
			ind = ind + 2**bit
			if (ind>n) then !Adding this bit makes ind>n, find the next smallest bit that fits
				ind = ind - 2**bit
				if (ind==n)return
				bit = int(log((n-ind)*1.)/log(2.))
				ind = ind + 2**bit
			end if
		else if (accum==tgt) then
			best = ind
			return
		else
			best = -1 !we somehow encountered an error
			return
		end if
	end do	
end function bitsearch

! Update the Fenwick tree due to a single update in the parent array.
! Use the delta between old/new array entry. Cost is O(log N).
pure subroutine update_fen(tree, delta, upd_index)
	integer, dimension(:), intent(inout) :: tree
	integer, intent(in) :: delta
	integer, intent(in) :: upd_index
	integer :: j
	
	j = upd_index
	do while (j<=size(tree))
		tree(j) = tree(j) + delta
		j = j + iand(j,-j)
	end do
end subroutine update_fen

! It is often the case that contiguous chunks of updates are made to the
! parent array. Naive updates would then have a O(M log N) cost for tree
! size N and M updates. Instead we can break these updates into two
! regions: inside the contiguous chunk and everything after. As a result
! we have reduced the cost to O(log N + M log M). For "large" M and N>>M
! this can translate to a significant savings, but even for M=2 or 3 it
! is still faster than 2 or 3 naive sequential updates.
pure subroutine update_fen_chunk(tree, delta, upd_index)
	integer, dimension(:), intent(inout) :: tree
	integer, dimension(:), intent(in) :: delta
	integer, intent(in) :: upd_index
	integer :: i,j,delta_tot,m

	m = size(delta)
	!update the tree part that occurs within the update chunk
	do i = 1, m
		j = upd_index + i-1
		do while (j<upd_index + m)
			tree(j) = tree(j) + delta(i)
			j = j + iand(j,-j)
		end do
	end do
	!update rest of tree past chunk, but now with total delta
	delta_tot = sum(delta)
	j = upd_index + m
	do while (j<=size(tree))
		tree(j) = tree(j) + delta_tot
		j = j + iand(j,-j)
	end do
end subroutine update_fen_chunk

! Compute the prefix sum for 1 to "ind" of the parent array. It is also
! possible to compute sums for index ranges subtracting the lower unused
! range from 1 to index start. Cost is O(log N).
pure function query_fen(tree, ind) result(psum)
	integer, dimension(:), intent(in) :: tree
	integer, intent(in) :: ind
	integer :: j, psum

	psum = 0
	j = ind
	do while (j>0)
		psum = psum + tree(j)
		j = j - iand(j,-j)
	end do
end function query_fen

end module fenwick_tree
