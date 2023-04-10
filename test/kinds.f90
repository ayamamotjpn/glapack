module kinds
	implicit none
#ifdef SINGLE
	integer,parameter :: DP=4
#else
	integer,parameter :: DP=8
#endif
	public DP
end module
