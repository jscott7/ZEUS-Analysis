c     ============================
      subroutine dt4(v1, v2, ans)
c     ============================
c Multiply two 4-vectors or the form (E - px - py - pz)
c
      implicit none
 
      real v1(4), v2(4)
      real ans

      ans =((v1(1)*v2(1))-(v1(2)*v2(2))-(v1(3)*v2(3))-(v1(4)*v2(4)))
 
      return
      end
