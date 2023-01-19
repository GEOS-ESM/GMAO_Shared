function driver (args)
field = subwrd  (args,1)
tag1  = subwrd  (args,2)
tag2  = subwrd  (args,3)
numf  = subwrd  (args,4)
title = subwrd  (args,5)

'set vpage off'
'set parea off'
'set mproj scaled'

* Mean Square Error
* -----------------
 position = '1 1 2 2'
'run statdplt 'field' Dmse 'tag1' 'tag2' 'numf' 'title' 'position

* Residual:  MSE - [Bias + Ampl + Phase]
* --------------------------------------
 position = '3 1 3 2'
'run statdplt 'field' Dres 'tag1' 'tag2' 'numf' 'title' 'position

* Bias Error
* ----------
 position = '1 2 3 2'
'run statdplt 'field' Dmes 'tag1' 'tag2' 'numf' 'title' 'position

* AMPL Error
* ----------
 position = '2 2 3 2'
'run statdplt 'field' Damp 'tag1' 'tag2' 'numf' 'title' 'position

* PHAZ Error
* ----------
 position = '3 2 3 2'
'run statdplt 'field' Dphz 'tag1' 'tag2' 'numf' 'title' 'position

'set vpage off'
'set parea off'

return
