pro test_fpraw_write,shotnum


;Th_lamp_468_calib_5m.nef
;restore,'/data/tatooine/FP/data/Calibration/0009114/Th_lamp_468_calib_5m_Calib.sav'

;multiple exposure test
;restore,'/data/tatooine/FP/data/0000001/0000001_savefile.sav'

restore,'/data/tatooine/FP/data/0010018/0010018_savefile.sav'

shotstr = string(shotnum,format='(I07)')

status = writefprawmdsplus(afp,shotnum)
print,status


end
