integer*4 flag_nodata,flag_lli,flag_good,flag_shrt,flag_no4,     &
          flag_lgjump,flag_lwjump,flag_gap,flag_nop2,flag_lowele,&
          flag_lwbad,flag_lgbad,flag_bigsd,flag_lccheck,         &
          flag_lwconn,flag_pcbad,flag_pc1ms

! c bit flag. (no flag means good data)
! 0 - 15 (ok), 16 -31 (del)
parameter(flag_nodata=31)
parameter(flag_no4   =30)
parameter(flag_lowele=29)
parameter(flag_shrt  =28)
parameter(flag_lwbad =25)
parameter(flag_lgbad =24)
parameter(flag_lccheck =23)
parameter(flag_pc1ms =18)
parameter(flag_pcbad =16)
parameter(flag_lwconn=6)
parameter(flag_lgjump=2 )
parameter(flag_lwjump=3 )
parameter(flag_gap   =4 )
parameter(flag_lli   =5 )
parameter(flag_bigsd =7 )
