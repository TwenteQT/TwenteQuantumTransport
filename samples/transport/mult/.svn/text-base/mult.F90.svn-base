Program mult
#define DEF_DBL_PREC 8
!!$Local
      Real (Kind=DEF_DBL_PREC) :: t_lr1 (2, 2) = 0.0d0, r_lr1 (2, 2) = 0.0d0, t_lr1_spec (2, 2) = 0.0d0
      Real (Kind=DEF_DBL_PREC) :: t_rl1 (2, 2) = 0.0d0, r_rl1 (2, 2) = 0.0d0, t_rl1_spec (2, 2) = 0.0d0
      Real (Kind=DEF_DBL_PREC) :: r_shc1 (2) = 0.0d0, l_shc1 (2) = 0.0d0
      Real (Kind=DEF_DBL_PREC) :: r_shc1x, l_shc1x
      Real (Kind=DEF_DBL_PREC) :: damp1r (2, 2) = 0.0d0, damp1i (2, 2) = 0.0d0
      Real (Kind=DEF_DBL_PREC) :: consv1, consw1
      Real (Kind=DEF_DBL_PREC) :: q1, q2, q3
      Real (Kind=DEF_DBL_PREC) :: p1, p2, p3
      Integer, Parameter :: fl = 200, fli = 201
      Character (Len=200) :: str


      Open (Unit=fli, File='andout', Action='read')


      Read (fli,*) str
      Read (fli, '("         ",2(5x,g13.6))') consv1, consw1
      Read (fli,*) str
      Read (fli, '("     by spin :",2(3x,f15.10))') l_shc1
      Read (fli, '("     total   :",3x,f15.10)') l_shc1x
      Read (fli,*) str
      Read (fli, '("     by spin :",2(3x,f15.10))') r_shc1
      Read (fli, '("     total   :",3x,f15.10)') r_shc1x
      Read (fli,*) str
      Read (fli,*) str
      Read (fli, '("   X    :",2(2x,g17.10))') damp1r (1, :)
      Read (fli, '("   Y    :",2(2x,g17.10))') damp1r (2, :)
      Read (fli,*) str
      Read (fli, '("   X    :",2(2x,g17.10))') damp1i (1, :)
      Read (fli, '("   Y    :",2(2x,g17.10))') damp1i (2, :)
      Read (fli,*) str
      Read (fli,*) str
      Read (fli,*) str
      Read (fli, '("   Up   :",2(3x,g17.10))') t_lr1 (1, :)
      Read (fli, '("   Down :",2(3x,g17.10))') t_lr1 (2, :)
      Read (fli,*) str
      Read (fli, '("         ",2(3x,g17.10))') p1, p2
      Read (fli, '("  Grand total:                  ",g17.10)') p3
      Read (fli,*) str
      Read (fli,*) str
      Read (fli, '("   Up   :",2(3x,g17.10))') t_rl1 (1, :)
      Read (fli, '("   Down :",2(3x,g17.10))') t_rl1 (2, :)
      Read (fli,*) str
      Read (fli, '("         ",2(3x,g17.10))') q1, q2
      Read (fli, '("  Grand total:                  ",g17.10)') q3
      Close (fli)

      l_shc1 = l_shc1 * 4.0d0
      r_shc1 = r_shc1 * 4.0d0
      l_shc1x = l_shc1x * 4.0d0
      r_shc1x = r_shc1x * 4.0d0
      t_lr1 = t_lr1 * 4.0d0
      t_rl1 = t_rl1 * 4.0d0
      q1=q1*4.0d0
      q2=q2*4.0d0
      q3=q3*4.0d0
      p1=p1*4.0d0
      p2=p2*4.0d0
      p3=p3*4.0d0
      
      Open (Unit=fl, File='andout_new', Action='write')

      Write (fl,*) ' Conservation: integrated,      worst'
      Write (fl, '("         ",2(5x,g13.6))') consv1, consw1

      Write (fl,*) ' Left Sharvin conductance :'
      Write (fl, '("     by spin :",2(3x,f15.10))') l_shc1
      Write (fl, '("     total   :",3x,f15.10)') l_shc1x
      Write (fl,*) ' Right Sharvin conductance :'
      Write (fl, '("     by spin :",2(3x,f15.10))') r_shc1
      Write (fl, '("     total   :",3x,f15.10)') r_shc1x
      Write (fl,*) '---------------------------------------------'
      Write (fl,*) ''
      Write (fl,*) ' Damping constants, Re part:'
      Write (fl, '("   X    :",2(2x,g17.10))') (damp1r(1, :))
      Write (fl, '("   Y    :",2(2x,g17.10))') (damp1r(2, :))
      Write (fl,*) ''
      Write (fl,*) ' Damping constants, Im part:'
      Write (fl, '("   X    :",2(2x,g17.10))') (damp1i(1, :))
      Write (fl, '("   Y    :",2(2x,g17.10))') (damp1i(2, :))
      Write (fl,*) '---------------------------------------------'
      Write (fl,*) ''

!!$       Write (fl,*) ' L->R Conductance(carefull with decomposition!):'
      Write (fl,*) ' L->R Conductance:'
      Write (fl,*) '              Up:                Down:'
      Write (fl, '("   Up   :",2(3x,g17.10))') t_lr1 (1, :)
      Write (fl, '("   Down :",2(3x,g17.10))') t_lr1 (2, :)
      Write (fl,*) ''
      Write (fl,*) ' Conductance total:'
      Write (fl, '("         ",2(3x,g17.10))') p1, p2
      Write (fl, '("  Grand total:                  ",g17.10)') p3
      Write (fl,*) ''
!!$       Write (fl,*) ' R->L Conductance(carefull with decomposition!):'
      Write (fl,*) ' R->L Conductance:'
      Write (fl,*) '              Up:                Down:'
      Write (fl, '("   Up   :",2(3x,g17.10))') t_rl1 (1, :)
      Write (fl, '("   Down :",2(3x,g17.10))') t_rl1 (2, :)
      Write (fl,*) ''
      Write (fl,*) ' Conductance total:'
      Write (fl, '("         ",2(3x,g17.10))') q1, q2
      Write (fl, '("  Grand total:                  ",g17.10)') q3
      Write (fl,*) ''
      Close (fl)
End Program
