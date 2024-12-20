module time_scheme

    use linear_algebra_mod
    implicit none

    public :: Advance

contains

    subroutine Advance(df, Un, tn, Unp1)

        !In
        type(DataType), intent(inout)         :: df
        real(pr), dimension(:), intent(in) :: Un
        real(pr), intent(in)               :: tn

        !Out
        real(pr), dimension(1:df%Nx*(df%jend-df%jbeg)), intent(out) :: Unp1

        !Local
        real(pr), dimension(1:df%Nx*(df%jend-df%jbeg)) :: Snp1

        Snp1 = SrcTermFunc(df, Un, tn+df%dt)
        !Snp1 = 1.0_pr
        
        call Lap_BiCGStab(df, Snp1, 10000, 1.e-8_pr, Unp1)

    end subroutine Advance


    

end module time_scheme