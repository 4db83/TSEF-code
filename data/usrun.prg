'caemp.prg
'pagedelete quarterly
wfcreate(wf=null) q 1947:1 2011:1

' data location
read(b14) .\GDPC1.xls x

' generate data transformations
genr t = @trend + 1
genr y = log(x)
genr dy = d(y)

smpl 1947:1 2011:1

' Enter the Non-ARMA terms. NOTE: This is a string of names that must be generated beforehand.
%nonarma = "dy c"

' Enter the maximum number of AR and MA terms you wish to use.
!ar = 1
!ma = 0

'Selecting models for y, within ARMA(4,4)
table(!ar+4,!ma+3) AIC
 AIC(2,3) = "MA Order"
 AIC(4,1) = "AR Order"
 AIC(1,1) = "AIC Values Various ARMA Models"

table(!ar+4,!ma+3) SIC
 SIC(2,3) = "MA Order"
 SIC(4,1) = "AR Order"
 SIC(1,1) = "SIC Values Various ARMA Models"

for !i=0 to !ar
	AIC(!i+4,2) = !i		
   	SIC(!i+4,2) = !i		
next
for !i=0 to !ma
	AIC(3,!i+3) = !i
   	SIC(3,!i+3) = !i
next   	

!f = !ar+1
!h = !ma+1

%ars = ""
%mas = ""

for !a =1 to !f
for !b =1 to !h

	for !j=1 to !a
		if !j -1 = 0 then
		%ar = ""
		else
		!q = !j-1
		%ar = %ars + " ar(" + @str(!q) +")"
		%ars = %ar
		endif
	next

	for !i=1 to !b
		if !i -1 = 0 then
		%ma = ""
		else
		!r = !i-1
		%ma = %mas + "ma(" + @str(!r) +")"
		%mas= %ma
		endif
	next

	!k = !a-1
	!m = !b-1
	equation arma{!k}{!m}.ls {%nonarma} {%ars} {%mas}

	AIC(!a + 3,!b + 2)=@aic
	SIC(!a + 3,!b + 2)=@schwarz
	
	%ars = ""
	%ar = ""
	%mas = ""
	%ma = ""
next
next


'show aic
'show sic

scalar psi  =  1/(1-(c(2)+c(3)+c(4)))
genr y_P    = psi*(y - c(2)*y(-1) - c(3)*y(-2) - c(4)*y(-3))
genr y_T    = y-y_P

y.hpf(lambda=1600) yP
genr yT = y-yP

group yy yT y_T
freeze(FIgure1) yy.line
show Figure1


stop

