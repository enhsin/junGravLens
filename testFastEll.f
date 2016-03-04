       program hello
          print *, "Hello World!"
c         implicit double precision (a-h, o-z)
c         double precision defl(2)
          call fastelldefl(0.1, 0.1, 0.2, 0.5, 0.6, 3, defl)	
          print *, defl(1), defl(2)
       end program hello