program main
    use gustavo
    use daniel
    
    real, dimension(5,5) :: A, B
    type(RowPacked) :: E
    type(ColPacked) :: C
    type(Pivot) :: Pivo
    
    A(1,:) = (/1.d0, 0.d0, 0.d0, -1.d0, 0.d0/)
    A(2,:) = (/2.d0, 0.d0, -2.d0, 0.d0, 3.d0/)
    A(3,:) = (/0.d0, -3.d0, 0.d0, 0.d0, 0.d0/)
    A(4,:) = (/0.d0, 4.d0, 0.d0, -4.d0, 0.d0/)
    A(5,:) = (/5.d0, 0.d0, -5.d0, 0.d0, 6.d0/)
    
    B(1,:) = (/1.d0, 0.d0, 2.d0, 3.d0, 0.d0/)
    B(2,:) = (/0.d0, 2.d0, 0.d0, 1.d0, 2.d0/)
    B(3,:) = (/2.d0, 0.d0, 1.d0, 0.d0, 0.d0/)
    B(4,:) = (/3.d0, 1.d0, 0.d0, 2.d0, 1.d0/)
    B(5,:) = (/0.d0, 2.d0, 0.d0, 1.d0, 1.d0/)
    
    E = GatherRow(B,MatrixDensity)
    
    
end program
