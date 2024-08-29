program greedyy
	!!!
	!
	! Main program of this file. Runs the tests - via the Test subroutine - then computes the fraction
	! that are successes, for each test, and writes those results to a specified file (named "GreedyOutputFile.txt").
	!
	!!!


    implicit none
    
    integer :: N, i, j, NoSuccess
    
    real, dimension(8,8) :: D
    
    real :: Frac
    
    real, dimension(19) :: results
    
    N = 8
    
    ! Loops running the specified tests, calculating the resulting fractions
    do i = 1, 19
        call Test(N, D, 0.05 * i, 100000, NoSuccess)
        results(i) = NoSuccess / 100000.0
    end do
    
	! Writes the results to the file    
        open (unit = 2, file = "GreedyOutputFile.txt")
        write (2,*) results
        close(2)
    
    end program greedyy
    
    
    subroutine Test(N, D, Frac, NoTests, NoSuccess)
	!!!
	!
	! Runs the tests for a single, specified, value of Frac
	!
	!!!
    
        integer, intent(in) :: N, NoTests
    
        real, intent(in) :: Frac
        
        real, intent(in), dimension(N, N) :: D	
    
        integer, dimension(N) :: OptimalTour, NearlyOptimalTour, GreedyResult
    
        integer :: NoSuccess
    
        real :: gap, MaxAlp, MaxErr, x
    
        logical :: Equi
    
    
        NoSuccess = 0
    
        do i = 1, NoTests
    
	    ! Randomly chooses parameters for generating an instance of TSP
            call random_number(x)
            MaxAlp = 6 + (10*x)
            call random_number(x)
            MaxErr = 3 + (MaxAlp - 4)*x
 
	    ! Calculates the required gap between the optimal and nearly optimal tours
            gap = ( 4 * Frac / (4 - Frac) ) * ( (MaxAlp + MaxErr) / 2 )
    
	    ! Generates a single instance of the TSP, with specified properties
            call ProbGen( N, MaxAlp, MaxErr, OptimalTour, NearlyOptimalTour, D, gap )
    
	    ! Solves (or attempts to the generated problem with a greedy algorithm)
            call greedy(D, N, GreedyResult)

	    ! Checks if the greedy algorithm has correctly solved the problem     
            if( Equi(OptimalTour, GreedyResult, N) ) then
                NoSuccess = NoSuccess + 1
            end if
    
        end do
    end subroutine Test
    
    
    
    subroutine greedy( D, N, GreedyResult )
	!!!
	!
	! An implementation of a greedy algorithm to heuristically solve instances of TSP
	!
	!!!

        implicit none
    
        integer, intent(in) :: N
    
        real, dimension(N, N) :: D
    
        integer, dimension(N) :: GreedyResult, pool
        
        integer :: curr1, curr2, LeftPos, RightPos, LeftChampPos, RightChampPos, i, j
    
        real :: current
    
        !!!!! Finds the shortest edge in the graph !!!!!
    
        !Sets initial values
    
        curr1 = 1
        curr2 = 2
    
        current = D(curr1, curr2)
    
        !Loops over all rows after the second row
        do i = 2, N
            !Loops over all elements in the selected row up to the diagonal
            do j = 1, (i-1)
                !If the element currently under consideration is smaller than any seen previously
                if ( D(i, j) < current ) then
    
                    curr1 = i
    
                    curr2 = j		
    
                    !Notes the value of the element
                    current = D(i, j)
            
                end if
            end do
        end do
    
        do i = 1, N
    
            pool(i) = i
    
        end do 
    
        call swap(pool, N, curr1, N)
        call swap(pool, N, curr2, N-1)
    
	! Once the shortest edge is found, the greedy solution is initialized containing only that edge
        GreedyResult(1) = curr1
        GreedyResult(2) = curr2
    
        LeftPos = 1
        RightPos = 2
    
	!!!!! Continually adds the shortest possible edge to form a complete tour !!!!!
    
        do i = 3,(N-1)
            
	    ! Finds the shortest edge that could be added to one end of the path (that will eventually be a tour) formed so far
            LeftChampPos = 1
            do j = 2, (N-i+1)
                if( D(pool(j), GreedyResult(LeftPos)) < D(pool(LeftChampPos), GreedyResult(LeftPos)))  then
                    LeftChampPos = j
                end if 
            end do
    
	    ! Finds the shortest edge that could be added to the other end of the path (that will eventually be a tour) formed so far
            RightChampPos = 1
            do j = 2, (N-i+1)
                if( D(pool(j), GreedyResult(RightPos)) < D(pool(RightChampPos), GreedyResult(RightPos))) then
                    RightChampPos = j
                end if 
            end do
    
	    ! Adds the shortest edge that could be added to either end of the path to the path
            if ( D(pool(LeftChampPos), GreedyResult(LeftPos)) < D(pool(RightChampPos), GreedyResult(RightPos)) )then
                GreedyResult( modulo(Leftpos - 2, N) + 1) = pool(LeftChampPos)
                call swap(pool, N, LeftChampPos, (N-i+1))
                LeftPos = modulo( LeftPos-2, N) + 1
            else
                GreedyResult( modulo(Rightpos, N) + 1) = pool(RightChampPos)
                call swap(pool, N, RightChampPos, (N-i+1))
                RightPos = modulo( RightPos, N) + 1
            end if
    
        end do
    
        GreedyResult(modulo(Rightpos, N) + 1) = pool(1)
    
    
    end subroutine greedy
    
    
    
    
    
    subroutine ProbGen( N, MaxAlpha, MaxError, OptimalTour, NearlyOptimalTour, D, gap )
	!!!
	!
	! Function to generate a single instance of the travelling salesman problem, with fine control over
	! which are the optimal and nearly optimal tours, and the difference in length between them.
	!
	!!!

        implicit none
    
        integer, intent(in)   :: N
    
        real, intent(in)      :: MaxAlpha, MaxError, gap
    
        integer, dimension(N) :: OptimalTour, NearlyOptimalTour
    
        real, dimension(N)    :: alphas
    
        real, dimension(N, N) :: D
    
        real 		      :: x
    
        integer               :: i, j, l, k, p, q
    
        
	! Loops over every 
        do i = 1, N
    
	    ! Randomly initializes the value of the alphas
            call random_number(x)
            alphas(i) = MaxAlpha * x
  
	    ! Initializes each edge as - temporarily - having zero length
            do j = 1, N
                D(i, j) = 0.0
            end do 
    
        end do
    
	! Randomly shuffles the optimal tour variable to choose what will be the optimal tour
        call FYShuffle(N, OptimalTour)
    
	! Sets the length of all edges in the optimal tour
        do k = 1, N
            l = modulo(k, N) + 1
            D(OptimalTour(k), OptimalTour(l)) = alphas(OptimalTour(k)) + alphas(OptimalTour(l))
            D(OptimalTour(l), OptimalTour(k)) = D(OptimalTour(k), OptimalTour(l))
        end do
    
	! Randomly chooses the nearly optimal tour
        call random_number(x)
        p = floor( N * x ) + 1
        q = modulo( P, N ) + 1
        NearlyOptimalTour = OptimalTour
        call swap(NearlyOptimalTour, N, p, q)
    
        ! Sets the length of the remaining edges with unassigned lengths in the nearly optimal tour
        D( OptimalTour(modulo(p-2, N) + 1), OptimalTour(q) ) = alphas(OptimalTour(modulo(p-2, N) + 1)) + alphas(q) + (gap/2)
        D( OptimalTour(q), OptimalTour(modulo(p-2, N) + 1) ) = D( OptimalTour(modulo(p-2, N) + 1), OptimalTour(q) )
        D( OptimalTour(modulo(q, N) + 1), OptimalTour(p) ) = alphas(OptimalTour(modulo(q, N) + 1)) + alphas(p) + (gap/2)
        D( OptimalTour(p), OptimalTour(modulo(q, N) + 1) ) = D( OptimalTour(modulo(q, N) + 1), OptimalTour(p) )
    
    
	! Randomly assigns lengths to the remaining edges in the generated TSP instance
        do i = 1, N
            do j = 1, N
                if( D(i, j) < 0.000000000001 ) then 
                    call random_number(x)
                    D(i, j) = alphas(i) + alphas(j) + (gap/2) + (MaxError - (gap/2))*x
                    D(j, i) = D(i, j)
                end if
            end do 
        end do
    
    end subroutine ProbGen
    
    
    
    
    subroutine swap(A, N, p, q)
	!!!
	!
	! A simple function to swap the position of any two elements of a given array
	!
	!!!

        implicit none
    
        integer :: l
        integer, intent(in) :: N, p, q
        integer, dimension(N) :: A
    
        !Selects the pth element of the given array, A, and stores it temporarily
        l = A(p)
        !Swaps the position of the pth and qth elements in the given array, A 
        A(p) = A(q)
        A(q) = l
    
    end subroutine
    
    
    
    SUBROUTINE FYShuffle(N, Shuffle)
    !!!
    !
    ! A function to generate a random tour of length, N, using a Fisher-Yates shuffle
    !
    !!!

    implicit none
    
    REAL :: X
    
    INTEGER :: i, l, N, p
    
    INTEGER, DIMENSION(N) :: Setup, Shuffle
    
    !Sets up an array of integers 1 to N
    DO i=1,N
        Shuffle(i) = i
    END DO
    
    
    !Do loop counting down from N to 1
    DO i=N,1,-1
    
        !Generates a random number, x, in (0, 1]
        CALL RANDOM_NUMBER(X)
        !Generates a random integer in (0,i) using x, to randomly select an integer not yet drawn
        p = FLOOR( (i)*x ) + 1
    
        call swap(Shuffle, N, p, i)
    
    END DO
    
    end subroutine FYShuffle
    
    logical function Equi(A, B, N) result(equiva)
	!!!
	!
	! A function to test if two tours are equivalent. Takes the two tours: A and B, plus their shared length as arguments
	!
	!!!
    
        integer :: i, j, m, l, counter 	
        
        integer, intent(in) :: N
    
        integer, dimension(N), intent(in) :: A, B	
    
        counter = 0
    
        !Counting the shared edges of the two tours
        !Loops over all nodes in tour A
        DO i=1,N
            !Loops over all nodes in tour B
            DO j=1,N
                !For each element in A finds the matching element in B
                IF( A(i) == B(j) ) THEN
                    !Calculates the next node in tour A
                    k = MODULO(i, N) + 1
                    !Calculates the next node in tour B
                    m = MODULO(j, N) + 1
                    !Calculates the previous node in tour B
                    l = MODULO(j-2, N) + 1
                    !If the next node in tour A is the same as the next or previous node in tour B
                    IF( A(k) == B(m) .or. A(k) == B(l) ) THEN
                        !Iterates the counter variable to mark that this edge is shared
                        counter = counter + 1
                    END IF
                END IF
            END DO
        END DO
    
        !Checks if all the edges in the two tour are shared 
        if ( counter == N ) then 
            equiva = .true.
        else
            equiva = .false.
        end if
    
    end function Equi
