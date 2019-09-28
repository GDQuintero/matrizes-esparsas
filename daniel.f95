module daniel

	use gustavo
    
    contains
    
!================================================================================================
!     PERMUTA DUAS COLUNAS DUMA MATRIZ EMPACOTADA POR **COLUNAS**, J1 É MENOR DO QUE J2
!================================================================================================
    subroutine col_perm_colpacked(A,j1,j2)! supposed that j1 < j2
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
		do i = 1, size(A%len_col)!pode-se melhorar este ciclo do, nao precisa percorrer o vetor len_col tudo
            A%col_start(i+1) = A%len_col(i)+A%col_start(i)
		enddo
    end subroutine col_perm_colpacked
!================================================================================================
!     PERMUTA DUAS LINHAS DUMA MATRIZ EMPACOTADA POR **COLUNAS**, I1 É MENOR DO QUE I2
!================================================================================================    
    subroutine row_perm_colpacked(A, i1, i2)
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
    end subroutine row_perm_colpacked
!================================================================================================
!     PERMUTA DUAS LINHAS DUMA MATRIZ EMPACOTADA POR **LINHAS**, I1 É MENOR DO QUE I2
!================================================================================================
    subroutine row_perm_rowpacked(A,i1,i2)! supposed that i1 < i2
        implicit none
        type(RowPacked) :: A
        integer :: i1, i2, j, m, i_aux
        integer, allocatable :: inter(:)
        real, allocatable :: array(:)
        real :: r_aux
        m = max(A%len_row(i1),A%len_row(i2))
        if(A%len_row(i1) /= A%len_row(i2)) then
            allocate(inter(m), array(m))
        end if
		if(A%len_row(i1) == A%len_row(i2)) then
			do j = 1, m
				r_aux = A%value(j-1+A%row_start(i2))
				i_aux = A%col_index(j-1+A%row_start(i2))
				A%value(j-1+A%row_start(i2)) = A%value(j-1+A%row_start(i1))
				A%col_index(j-1+A%row_start(i2)) = A%col_index(j-1+A%row_start(i1))
				A%value(j-1+A%row_start(i1)) = r_aux
				A%col_index(j-1+A%row_start(i1)) = i_aux
			enddo
		elseif(A%len_row(i1) > A%len_row(i2)) then
			do j = 1, m
				array(j) = A%value(j-1+A%row_start(i1))
				inter(j) = A%col_index(j-1+A%row_start(i1))
			enddo
			do j = 1, A%len_row(i2)!escrevemos as entradas da linha i2 no lugar da i1
                A%value(j-1+A%row_start(i1)) = A%value(j-1+A%row_start(i2))
                A%col_index(j-1+A%row_start(i1)) = A%col_index(j-1+A%row_start(i2))
			enddo
			if(i2-i1 > 1) then!empurramos à esquerda os elementos entre as linhas i1 e i2
                do j = 1, A%row_start(i2)-A%row_start(i1+1)
                    A%value(j-1+A%row_start(i1)+A%len_row(i2)) = A%value(j-1+A%row_start(i1+1))
                    A%col_index(j-1+A%row_start(i1)+A%len_row(i2)) = A%col_index(j-1+A%row_start(i1+1))
                enddo
            endif
            do j = 1, m!escrevemos as entradas da linha i1 no lugar da i2
                A%value(j-1+A%row_start(i1)+A%len_row(i2)+A%row_start(i2)-A%row_start(i1+1)) = array(j)
                A%col_index(j-1+A%row_start(i1)+A%len_row(i2)+A%row_start(i2)-A%row_start(i1+1)) = inter(j)
            enddo
		else
			do j = 1, m
				array(j) = A%value(j-1+A%row_start(i2))
				inter(j) = A%col_index(j-1+A%row_start(i2))
			enddo
			if(i2 == size(A%len_row)) then 
                do j = A%len_row(i1), 1, -1
                    A%value(j+size(A%value)-A%len_row(i1)) = A%value(j-1+A%row_start(i1))
                    A%col_index(j+size(A%value)-A%len_row(i1)) = A%col_index(j-1+A%row_start(i1))
                enddo
            else !da para colocar elseif(i2 < m) mas isso é evidente
                do j = A%len_row(i1), 1, -1
                    A%value(j-A%len_row(i1)+A%row_start(i2+1)-1) = A%value(j-1+A%row_start(i1))
                    A%col_index(j-A%len_row(i1)+A%row_start(i2+1)-1) = A%col_index(j-1+A%row_start(i1))
                enddo
			endif
			if(i2-i1 > 1) then!empurramos de novo
                do j = A%row_start(i2)-A%row_start(i1+1), 1, -1
                    A%value(j+A%row_start(i1+1)-1+A%len_row(i2)-A%len_row(i1)) = A%value(j+A%row_start(i1+1)-1)
                    A%col_index(j+A%row_start(i1+1)-1+A%len_row(i2)-A%len_row(i1)) = A%col_index(j+A%row_start(i1+1)-1)
                enddo
			endif
			do j = 1, m
                A%value(j-1+A%row_start(i1)) = array(j)
                A%col_index(j-1+A%row_start(i1)) = inter(j)
			enddo
		endif
		i_aux = A%len_row(i2)
		A%len_row(i2) = A%len_row(i1)
		A%len_row(i1) = i_aux
		do j = 1, size(A%len_row)!pode-se melhorar este ciclo do, nao precisa percorrer o vetor len_row tudo
            A%row_start(j+1) = A%len_row(j)+A%row_start(j)
		enddo
    end subroutine row_perm_rowpacked
!================================================================================================
!     PERMUTA DUAS COLUNAS DUMA MATRIZ EMPACOTADA POR **LINHAS**, j1 É MENOR DO QUE j2
!================================================================================================    
    subroutine col_perm_rowpacked(A, j1, j2)
        implicit none
        type(RowPacked) :: A
        integer :: j1, j2, i
        do i = 1, size(A%value)
            if(A%col_index(i) == j1) then
                A%col_index(i) = j2
            elseif(A%col_index(i) == j2) then
                A%col_index(i) = j1
            endif
        enddo
    end subroutine col_perm_rowpacked
!================================================================================================
!     ENCONTRA O PIVOTE OTIMO SEGUNDO O CRITERIO DE MARKOWITZ, MATRIZ EMPACOTADA POR **COLUNAS**
!================================================================================================  
	function Markowitz(A)
		implicit none
		type(pivot) :: Markowitz
		type(ColPacked) :: A
		integer, allocatable :: aux(:,:), lrow(:)!lrow es un vector auxiliar interno
		real, parameter :: u = 0.25!threshold pivoting
		integer :: i, j, k, l, n, tau, min_prod
		logical :: bool
		real :: r
		n = size(A%len_col)!numero das colunas de A
		tau = size(A%row_index)
		allocate(aux(n,n), lrow(n))
		lrow = 0
		do i = 1, tau
            lrow(A%row_index(i)) = lrow(A%row_index(i)) + 1
		enddo
		do i = 1, n
            aux(i,:) = (lrow(i)-1)*(A%len_col-1)
		enddo
		min_prod = minval(aux)
		if(min_prod == 0) then
            do i = 1, n
                if(A%len_col(i) == 1) then
                    Markowitz%col = i
                    Markowitz%row = A%row_index(A%col_start(i))
                    return
                elseif(lrow(i) == 1) then
                    do j = 1, n
                        do k = A%col_start(j), A%col_start(j+1) 
                            if(A%row_index(k) == i) then
                                Markowitz%col = j
                                Markowitz%row = i
                                return
                            endif
                        enddo
                    enddo
                endif
            enddo
        else
            do j = 1, n
                do i = 1, n
                    if(aux(i,j) == min_prod) then
                        k = A%col_start(j)
                        bool = .true.
                        do while(bool .and. k < A%col_start(j+1))
                            if(A%row_index(k) == i) then 
                                r = abs(A%value(k))
                                bool = .false.
                            endif
                            k = k + 1
                        enddo
                        l = A%col_start(j)
                        bool = .true.
                        do while(bool .and. l < A%col_start(j+1) .and. l /= (k-1))!threshold pivoting
                            if(r < u*abs(A%value(l))) then
                                bool = .false.
                                exit!dale una revisada daniel del furuto
                            endif
                            l = l + 1
                        enddo
                        if(bool) then
                            Markowitz%col = j
                            Markowitz%row = i
                            return!que acontece cuando no consigue ni un pivo numericamente estable? que valor retorna?
                        endif
                    endif
                enddo
            enddo
		endif
	end function
!================================================================================================
!     SIMPLER STRATEGIES MIN ROW IN MIN COL PARA COLPACKED
!================================================================================================  
	function MinRowInMinCol(A)
		type(pivot) :: MinRowInMinCol
		type(ColPacked) :: A
		integer :: i, j, k, l, m, n, min_row, min_col, tau
		integer, allocatable :: len_row(:)
		logical :: bool
		real, parameter :: u = 0.25!threshold pivoting
		n = size(A%len_col)
		tau = size(A%row_index)
		allocate(len_row(n))
		len_row = 0!tal vez innecesario
		do i = 1, tau ! creamos el vector len_row
            len_row(A%row_index(i)) = len_row(A%row_index(i)) + 1
		enddo
		min_col = minval(A%len_col)
		do j = 1, n
			if(A%len_col(j) == min_col) then
				if(len_row(A%row_index(A%col_start(j))) .le. len_row(A%row_index(A%col_start(j)+1))) then!las primeras 2 entradas de la col j
					i = A%row_index(A%col_start(j))
				else
					i = A%row_index(A%col_start(j)+1)
				endif
				do k = A%col_start(j)+2, A%col_start(j+1)-1
					if(len_row(A%row_index(k)) < len_row(i)) then! .and. (abs(A%value(k)) .ge. u*abs(A%value(k+1)))) then
						i = A%row_index(k)
					endif
				enddo
				l = A%col_start(j)
				bool = .true.
				do while(bool .and. l < A%col_start(j+1))!threshold pivoting
					if(abs(A%value(i)) < u*abs(A%value(l))) then
						bool = .false.
						exit!dale una revisada daniel del furuto
					endif
					l = l + 1
				enddo
				if(bool) then
					MinRowInMinCol%col = j	
					MinRowInMinCol%row = i
					MinRowInMinCol%value = A%value(i)
					return
				endif
			endif
		enddo
	end function
!================================================================================================
!     IMPRIMIENDO MATRIZ EMPAQUETADA
!================================================================================================  
    function printf(A)
        implicit none
        type(colpacked) :: A
        character(len = 3), allocatable :: printf(:,:)
        integer :: i, j, n
        n = size(A%len_col)
        allocate(printf(n,n))
        printf = "   "
        do j = 1, n
            do i = A%col_start(j), A%col_start(j+1)-1
                if(abs(A%value(i)) .ge. 10e-8) then
                    printf(A%row_index(i),j) = " x "
                else
                    printf(A%row_index(i),j) = "   "
                endif
            enddo
        enddo
        
        do i = 1, n
            print*, "|", printf(i,:), "|"
        enddo
	end function
end module
