SUBROUTINE my_sum_sub(x,y,out)
  REAL,intent(out):: out
  REAL,intent(in):: x, y
  out = x + y
  WRITE(*,*) 'x = ', x
  WRITE(*,*) 'y = ', y
  WRITE(*,*) 'sum = ', out
END SUBROUTINE
