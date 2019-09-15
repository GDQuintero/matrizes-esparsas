module daniel

	use gustavo
    
    contains
    
!================================================================================================
!     PERMUTA DUAS COLUNAS DUMA MATRIZ EMPACOTADA POR COLUNAS, J1 É MENOR DO QUE J2
!================================================================================================
    subroutine col_permutation(A,j1,j2)! supposed that j1 < j2
        implicit none
        type(ColPacked) :: A
        integer :: j1, j2, i, m, i_aux
        integer, allocatable :: inter(:)
        real, allocatable :: array(:)
        real :: r_aux
        m = max(A%len_col(j1),A%len_col(j2))
        if(A%len_col(j1) /= A%len_col(j2)) then
            allocate(inter(m), array(m))
        end if
		if(A%len_col(j1) == A%len_col(j2)) then
			do i = 1, m
				r_aux = A%value(i-1+A%col_start(j2))
				i_aux = A%row_index(i-1+A%col_start(j2))
				A%value(i-1+A%col_start(j2)) = A%value(i-1+A%col_start(j1))
				A%row_index(i-1+A%col_start(j2)) = A%row_index(i-1+A%col_start(j1))
				A%value(i-1+A%col_start(j1)) = r_aux
				A%row_index(i-1+A%col_start(j1)) = i_aux
			enddo
		elseif(A%len_col(j1) > A%len_col(j2)) then
			do i = 1, m
				array(i) = A%value(i-1+A%col_start(j1))
				inter(i) = A%row_index(i-1+A%col_start(j1))
			enddo
			do i = 1, A%len_col(j2)!escrevemos as entradas da coluna j2 no lugar da j1
                A%value(i-1+A%col_start(j1)) = A%value(i-1+A%col_start(j2))
                A%row_index(i-1+A%col_start(j1)) = A%row_index(i-1+A%col_start(j2))
			enddo
			if(j2-j1 > 1) then!empurramos à esquerda os elementos entre as colunas j1 e j2
                do i = 1, A%col_start(j2)-A%col_start(j1+1)
                    A%value(i-1+A%col_start(j1)+A%len_col(j2)) = A%value(i-1+A%col_start(j1+1))
                    A%row_index(i-1+A%col_start(j1)+A%len_col(j2)) = A%row_index(i-1+A%col_start(j1+1))
                enddo
            endif
            do i = 1, m!escrevemos as entradas da coluna j1 no lugar da j2
                A%value(i-1+A%col_start(j1)+A%len_col(j2)+A%col_start(j2)-A%col_start(j1+1)) = array(i)
                A%row_index(i-1+A%col_start(j1)+A%len_col(j2)+A%col_start(j2)-A%col_start(j1+1)) = inter(i)
            enddo
		else
			do i = 1, m
				array(i) = A%value(i-1+A%col_start(j2))
				inter(i) = A%row_index(i-1+A%col_start(j2))
			enddo
			if(j2 == size(A%len_col)) then 
                do i = A%len_col(j1), 1, -1
                    A%value(i+size(A%value)-A%len_col(j1)) = A%value(i-1+A%col_start(j1))
                    A%row_index(i+size(A%value)-A%len_col(j1)) = A%row_index(i-1+A%col_start(j1))
                enddo
            else !da para colocar elseif(j2 < m) mas isso é evidente
                do i = A%len_col(j1), 1, -1
                    A%value(i-A%len_col(j1)+A%col_start(j2+1)-1) = A%value(i-1+A%col_start(j1))
                    A%row_index(i-A%len_col(j1)+A%col_start(j2+1)-1) = A%row_index(i-1+A%col_start(j1))
                enddo
			endif
			if(j2-j1 > 1) then!empurramos de novo
                do i = A%col_start(j2)-A%col_start(j1+1), 1, -1
                    A%value(i+A%col_start(j1+1)-1+A%len_col(j2)-A%len_col(j1)) = A%value(i+A%col_start(j1+1)-1)
                    A%row_index(i+A%col_start(j1+1)-1+A%len_col(j2)-A%len_col(j1)) = A%row_index(i+A%col_start(j1+1)-1)
                enddo
			endif
			do i = 1, m
                A%value(i-1+A%col_start(j1)) = array(i)
                A%row_index(i-1+A%col_start(j1)) = inter(i)
			enddo
		endif
		
		i_aux = A%len_col(j2)
		A%len_col(j2) = A%len_col(j1)
		A%len_col(j1) = i_aux
		
		do i = 1, size(A%len_col)
            A%col_start(i+1) = A%len_col(i)+A%col_start(i)
		enddo
    end subroutine col_permutation
!================================================================================================
!     PERMUTA DUAS LINHAS DUMA MATRIZ EMPACOTADA POR COLUNAS, I1 É MENOR DO QUE I2
!================================================================================================    
    subroutine row_permutation(A, i1, i2)
        implicit none
        type(ColPacked) :: A
        integer :: i1, i2, i
        do i = 1, size(A%value)
            if(A%row_index(i) == i1) then
                A%row_index(i) = i2
            elseif(A%row_index(i) == i2) then
                A%row_index(i) = i1
            endif
        enddo
    end subroutine

end module
