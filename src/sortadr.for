	program sortadr
c $Id: sortadr.for 181 2008-08-14 23:19:51Z tjansson $
c Program sortadr programmed by C.C.Tscherning, Geophysical
c Institute University of Copenhagen, Copyright 1992.
comment Program sort_adr, udarbejdet af C.C.Tscherning,
cfeb. 1979. GI reg. no. 82110. Last change 1995.03.21 by cct.
c
c Instructions in English and Danish:
c The program will sort items from a list (file) using a list of
c identification numbers. The file may contain anything like adresses
c or litterature references. An item may have several identification
c numbers, such as a unique id. and a group code.
C
C Each item must start with a  < followed by a text-string used for
c sorting purposes, followed by the id. numbers, terminated by -1.
C Example:
c       <tscherning  10055  3.90 4.111 -1
c The whole list is terminated by giving a start id. with just
c         <   -1
C
c Programmet sorterer adresser ud fra en adresseliste ved
c hj{lp af en sorteringsliste indeholdende identifika-
c tionsnumre p} adresserne. Listen kan selvf|lgelig godt inde-
c holde andet end adresser, f.eks. litteratur referancer.
c
c En adresse kan have flere numre. F.eks. et l|benummer
c samt et tal der angiver at adressen er p} en person eller
c institution, der tilh|rer en bestemt gruppe, (en IAG special
c studiegruppe f.eks.).
c
c Listen af identifikationsnumre skal begynde med en tekststring
c hvis f|rste element er et < (f.eks. <tscherningcc ), og skal
c terminere med et negativt tal (-1). Listen m} maximalt inde-
c holde 20 elementer. Herefter f|lger s} adressen, der maximalt
c m} v{re fordelt over 6 linier.
c
c Eksempel:
c _     <tscherningcc 22 5.26  3.30 -1
c _     C.C.Tscherning,
c _     etc. -
c
c Hele listen af adresser termineres ved hj{lp af en 'tom'
c liste af identifikations-numre, dvs.  <  -1 .
c
c Den indledende tekststring ( <text ) kan benyttes, hvis listen
c |nskes ordnet alfabetisk, se programmet sortlisttx, GI reg. no.
c 82111.
c
c Ved sorteringen benyttes en liste af numre, i det fgl. kaldet
c sorteringslisten. Denne indeholder identifikationsnumre
c p} de adresser, der skal udskrives. Den termineres af et
c negativt tal og m} maximalt indeholde 500 numre.
c
c \nskes hele listen udskrevet angives en tom sorteringsliste,
c dvs. blot best}ende af et negativt tal.
c
c Adresserne skal st} i et dokument p} baggrungslageret. Ved
c kaldet af dette program skal dokumentnavnet angives f|rst.
c
c Eksempler:
c Programkald:    sortadr
c Adresse dok.:   cctadr
c Papir type  :   0
c Sort.liste:     5.26 102  103 110 160 -1
c giver adresser p} medlemmer af studiegruppe 5.26 og p}
c 4 enkeltpersoner/institutioner.
c Programkald:    sortadr
c Adresse dok.:   song79adr
c Sort.liste  :   -1
c
c Bem{rk at indholdet af identifikations eller sorteringslisterne
c ikke beh|ver at v{re ordnet.
comment indl{sning af navnet p} adressefilen og papir type
        implicit real*8 (a-h,o-z)
	logical cont,EOL,EOF,wpout
        real sortlist(500),num(150)
	dimension iline(80)
	character*72 resfile,input
        CHARACTER*1 ITEM(80),RECO(512)
        CHARACTER*80 LINE
        EQUIVALENCE (LINE,ITEM(1))
        EOL=.TRUE.
	wpout=.FALSE.
	write(6,*)' Program Sortadr, Ver. feb. 92.'
	write(6,*)' Input name of inputfile '
	read(5,5) input
   5    format(a12)
	open(10,FILE=input)
	write(6,*)' Input name of output file '
	read(5,5) resfile
	open(20,FILE=resfile)
c ls sorteringsliste
	write(6,*)
     *' Input list of item numbers terminating with neg. nu.'
	m= 0
   10   continue
	ai=readno(5,EOL,EOF)
c       write(6,999)AI
c 999   format(D16.1)
        wpout= (abs(ai+2.0).lt.1.0E-8)
	if (wpout) write(*,*)' output in wp mode '
	if (ai.lt.0.0) goto 20
	m= m+1
        if (m.gt.500) stop ' m too large '
	sortlist(m)= ai
	goto 10
   20   continue
        kka=0
c l{s sorteringskoder
c begin j loop
        j=0 
        n=1 
  30    if (j.gt.0) backspace 10
    	AJ=readNO(10,EOL,EOF)
        j=j+1 
   	if (aj.lt.0.0) goto 80
	n= 1
	num(1)= aj
c begin k loop
   40   AJ=readno(10,EOL,EOF)
	if (aj.lt.0.0) goto 50
	n= n+1
c 	write(*,*)n,aj  
	num(n)= aj
        goto 40
c end k loop
  50    continue
c
	cont=m.gt.0
	i= 0
c BEGIN i loop
 55     I=I+1
	if (cont.AND.I.LE.M) THEN
          K=0
	  nchar=0
c BEGIN k loop
 60       K=K+1
          IF (CONT.AND.K.LE.N) THEN
            cont=abs(num(k)-sortlist(i)).gt.0.0e-6
            GOTO 60
          ENDIF
          GOTO 55
c end k loop
        ENDIF
c end i loop
c
c k loop, l{s adressen linievis
  70   continue
	read(10,'(A)')line
	ilen=0
c       read(10,72)line
  72    format(a72)
	if (ITEM(1).eq.'<'.OR.ITEM(2).EQ.'<'.OR.ITEM(3).EQ.'<') THEN
	if ((.not.cont).and.wpout)
     *  write(20,773)(RECO(kk),kk=1,nchar)
        if (.not.cont) kka=kka+1
	IF ((.not.cont).and.wpout)write(6,*)nchar
	nchar=0 
	goto 30
	end if
	if (cont) go to 70 
	do 772 kk=1,80
         iline(kk)=ichar(line(kk:kk))
 772    if (iline(kk).ne.32) iend=kk
	nchar1=nchar+1
	if (nchar.gt.0)RECO(nchar1)=' '
	do 774 kk=1,iend
	kk1=kk+nchar1
 774    RECO(kk1)=LINE(KK:KK)
	nchar=nchar1+iend
 773    format(256A1)
         if (.not.cont.and.(.not.wpout)) then
        write(20,773)(LINE(KK:KK),KK=1,iend)
        end if
        GOTO 70
c
c end j loop
  80    continue
c
	close(10)
c 26=end of medium
	i=26
	write(20,901)i
  901  format(a2)
     	close(20)
        write(6,902)j
  902   format(' items read =',i6) 
        write(6,903)kka
  903   format(' items written ',I6)
	stop
	end
C
C
      INTEGER FUNCTION IREADCHAR(LU,C,EOL,EOF)
c ireadchar return             the line position
c lu       call               the logical input unit, opened as sequential
c                             and formatted
c c        return             the next character from "lu"
c eol      call and return    at call: true, if next line is wanted
c                             at return: true, if lpos=last character in line
c eof     return              true, if end of file
c note: "eol" must be true at first call
      LOGICAL EOL,EOF
      CHARACTER*1 C,SPACE
      CHARACTER*1 ITEM(132)
      CHARACTER*132 LINE
      COMMON LAST,LPOS,LINE
      EQUIVALENCE (LINE,ITEM(1))
      SPACE=' '
      EOF=.FALSE.
  10  IF (EOL) THEN
        LPOS=0
        EOL=.FALSE.
        READ(LU,200,END=300)LINE
c       write(6,*)'linie ',(item(k),k=1,10)
c       write(6,*)LINE
c       write(6,999)(ITEM(IGG),IGG=1,6)
c999    format(6a1)
        LAST=133
        DO 12 LST=1,132
         LAST=LAST-1
         IF (ITEM(LAST).NE.SPACE) GOTO 14
  12    CONTINUE
  14    CONTINUE
      ENDIF
      LPOS=LPOS+1
      C=ITEM(LPOS)
      IREADCHAR=LPOS
      IF (LPOS.GE.LAST) EOL=.TRUE.
 200  FORMAT(A132)
      GOTO 400
 300  EOF=.TRUE.
 400  CONTINUE
      RETURN
      END
C
C
      DOUBLE PRECISION FUNCTION READNO(LU,EOL,EOF)
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL EOL,EOF
C function   read no
c -------------------
c read a numeric variable as text from the logical unit lu
c and convert to double precision real (preceeding text is skipped)
c at first call eol must be true
c a no must be first in line or preceeded by a space
c the function uses the "ireadchar" function
      CHARACTER*1 CHAR,SPACE,LASTCHAR
      LOGICAL NEG,FIRST,LEGAL
      SPACE=' '
  1   FIRST=.TRUE.
      IDEC=0
      NEG=.FALSE.
      R=0.0D0
  10  IF (EOF) GOTO 99
      IF (EOL) CHAR=SPACE
      LASTCHAR=CHAR
      III=IREADCHAR(LU,CHAR,EOL,EOF)
c     write(6,111)iii,char,lastchar
c 111  format(' readno ',a1,1x,a1,1x,a1)
      LEGAL=(CHAR.GE.'0'.AND.CHAR.LE.'9')
     .      .OR.CHAR.EQ.'.'.OR.CHAR.EQ.'+'.OR.CHAR.EQ.'-'
      IF (LEGAL.AND.LASTCHAR.EQ.SPACE) GOTO 20
      GOTO 10
   20 CONTINUE
c read and save sign, if any
      NEG = CHAR.EQ.'-'
      IF (CHAR.EQ.'-'.OR.CHAR.EQ.'+') THEN
        IF (EOL.OR.EOF) GOTO 1
        III=IREADCHAR(LU,CHAR,EOL,EOF)
      ENDIF
c digits before decimal point
   30 IF (CHAR.GE.'0'.AND.CHAR.LE.'9') THEN
        FIRST=.FALSE.
        I=ICHAR(CHAR)
        I=I - 48
        R=1.0D0*I+R*10.0D0
        IF (EOL.OR.EOF) GOTO 90
        III=IREADCHAR(LU,CHAR,EOL,EOF)
        GOTO 30
      ENDIF
c if digits after decimal point
      IF (CHAR.EQ.'.') THEN
 40     IF (EOL.OR.EOF) GOTO 90
        III=IREADCHAR(LU,CHAR,EOL,EOF)
        IF (CHAR.GE.'0'.AND.CHAR.LE.'9') THEN
          FIRST=.FALSE.
          I=ICHAR(CHAR)
          I=I - 48
          R= R*10.0D0 + I
          IDEC=IDEC+1
          GOTO 40
        ENDIF
      ENDIF
  90  IF (FIRST) GOTO 1
      DO 50 I=IDEC,1,-1
      R=R/10.0D0
   50 CONTINUE
      IF (NEG) R=-R
c      write(6,*)' r = ',r
   99 READNO=R
      RETURN
      END
