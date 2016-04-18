FUNCTION my_sum(x,y) result(out)
  REAL:: out
  REAL:: x, y
  WRITE(*,*) 'x = ', x
  WRITE(*,*) 'y = ', y
  out = x + y
END
