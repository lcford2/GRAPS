c  THIS SOFTWARE MAY NOT BE COPIED TO MACHINES OUTSIDE THE SITE FOR
c  WHICH IT HAD BEEN PROVIDED.  SEE "Conditions for External Use"
c  BELOW FOR MORE DETAILS.  INDIVIDUALS INTERESTED IN OBTAINING
c  THE SOFTWARE SHOULD CONTACT THE AUTHORS.
c
      subroutine FFSQP(nparam,nf,nineqn,nineq,neqn,neq,mode,iprint,
     *                  miter,inform,bigbnd,eps,epseqn,udelta,bl,bu,x,
     *                  f,g,iw,iwsize,w,nwsize,obj,constr,gradob,gradcn)
c
c     implicit real*8(a-h,o-z)
      integer nparam,nf,neqn,nineqn,nineq,neq,mode,iprint,miter,inform,
     *        iwsize,nwsize
      integer iw(iwsize)
      double  precision bl(1),bu(1),x(1),
     *        f(1),g(1),w(1)
c     double  precision bl(nparam),bu(nparam),x(nparam),
c    *        f(nf),g(nineq+neq),w(nwsize)
c
c   When nf=0 and/or nineq+neq=0, the dimension for f and/or g
c   should be at least 1 due to F77 standard. See more details
c   below.
c
      double  precision bigbnd,eps,epseqn,udelta
      external obj,constr,gradob,gradcn
c
c**********************************************************************c
c                                                                      c
c brief specification of various arrays and parameters in the calling  c
c sequence. See manual for more detailed description.                  c
c                                                                      c
c nparam : number of variables                                         c
c nf     : number of objective functions                               c
c nineqn : number of nonlinear inequality constraints                  c
c nineq  : number of inequality constraints                            c
c neqn   : number of nonlinear equality constraints                    c
c neq    : number of equality constraints                              c
c mode   : mode=CBA specifies job options as described below:          c
c          A = 0 : ordinary minimax problems                           c
c            = 1 : ordinary minimax problems with each individual      c
c                  function replaced by its absolute value, ie,        c
c                  an L_infty problem                                  c
c          B = 0 : monotone decrease of objective function             c
c                  after each iteration                                c
c            = 1 : monotone decrease of objective function after       c
c                  at most four iterations                             c
c          C = 1 : during line search, the function that rejected      c
c                  the previous step size is checked first;            c
c                  all functions of the same type ("objective" or      c
c                  "constraints") as the latter will then be checked   c
c                  first                                               c
c          C = 2 : all contraints will be checked first at every trial c
c                  point during the line search                        c
c iprint : print level indicator with the following options            c
c          iprint=0: no normal output except error information         c
c                    (this option is imposed during phase 1)           c
c          iprint=1:  a final printout at a local solution             c
c          iprint=2:  a brief printout at the end of each iteration    c
c          iprint=3:  detailed infomation is printed out at the end    c
c                     of each iteration for debugging purpose          c
c          iprint=10*N+M: N any positive integer, M=2 or 3.            c
c                     Information corresponding to iprint=M will be    c
c                     displayed at every 10*Nth iterations at the last c
c                     iteration                                        c
c miter  : maximum number of iterations allowed by the user to solve   c
c          the problem                                                 c
c inform : status report at the end of execution                       c
c          inform= 0:normal termination                                c
c          inform= 1:no feasible point found for linear constraints    c
c          inform= 2:no feasible point found for nonlinear constraints c
c          inform= 3:no solution has been found within miter iterations
c          inform= 4:stepsize is smaller than machine precision before c
c                    a successful new iterate is found                 c
c          inform= 5:failure of the QP solver in attempting to         c
c                    construct d0. A more robust QP solver may succeed.c
c          inform= 6:failure of the QP solver in attempting to         c
c                    construct d1. A more robust QP solver may succeed.c
c          inform= 7:inconsistent input data                           c
c          inform= 8:new iterate essentially identical to previous     c
c                    iterate, though stopping criteria not satisfied   c
c          inform= 9:penalty parameter too large, unable to satisfy    c
c                    nonlinear equality constraint                     c
c bigbnd : plus infinity                                               c
c eps    : stopping criterion that ensures at a solution, the norm of  c
c          the Newton direction vector is smaller than eps             c
c epseqn : tolerance of the violation of nonlinear equality constraintsc
c          allowed by the user at an optimal solution                  c
c udelta : perturbation size in computing gradients by finite          c
c          difference and the true perturbation is determined by       c
c          sign(x_i) X max{udelta, rteps X max{1, |x_i|}} for each     c
c          component of x, where rteps is the square root of machine   c
c          precision                                                   c
c bl     : array of dimension nparam,containing lower bound of x       c
c bu     : array of dimension nparam,containing upper bound of x       c
c x      : array of dimension nparam,containing initial guess in input c
c          and final iterate at the end of execution                   c
c f      : array of dimension max{1,nf}, containing objective values   c
c          at x in output                                              c
c g      : array of dimension max{1,nineq+neq}, containing constraint  c
c          values at x in output                                       c
c iw     : integer working space of dimension iwsize                   c
c iwsize : length of integer array iw                                  c
c w      : double precision working space of dimension nwsize.         c
c          at output, it contains lagrange multipliers                 c
c nwsize : length of double precision array w                          c
c obj    : subroutine that returns the value of objective functions    c
c          one upon each call                                          c
c constr : subroutine that returns the value of constraints            c
c          one upon each call                                          c
c gradob : subroutine that computes gradients of f, alternatively      c
c          it can be replaced by grobfd that computes finite           c
c          difference approximations                                   c
c gradcn : subroutine that computes gradients of g, alternatively      c
c          it can be replaced by grcnfd that computes finite           c
c          difference approximations                                   c
c                                                                      c
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
c                                                                      c
c                                                                      c
c                       FFSQP  Version 3.7b                            c
c                                                                      c
c                   Jian L. Zhou, Andre L. Tits,                       c
c                        and C.T. Lawrence                             c
c                   Institute for Systems Research                     c
c                               and                                    c
c                Electrical Engineering Department                     c
c                     University of Maryland                           c
c                     College Park, Md 20742                           c
c                                                                      c
c                         January, 1998                                c
c                                                                      c
c                                                                      c
c  The purpose of FFSQP is to solve general nonlinear constrained      c
c  minimax optimization problems of the form                           c
c                                                                      c
c   (A=0 in mode)     minimize    max_i f_i(x)   for i=1,...,n_f       c
c                        or                                            c
c   (A=1 in mode)     minimize    max_j |f_i(x)|   for i=1,...,n_f     c
c                       s.t.      bl   <= x <=  bu                     c
c                                 g_j(x) <= 0,   for j=1,...,nineqn    c
c                                 A_1 x - B_1 <= 0                     c
c                                                                      c
c                                 h_i(x)  = 0,   for i=1,...,neqn      c
c                                 A_2 x - B_2  = 0                     c
c                                                                      c
c                                                                      c
c                                                                      c
c                  Conditions for External Use                         c
c                  ===========================                         c
c                                                                      c
c   1. The FFSQP routines may not be distributed to third parties.     c
c      Interested parties should contact the authors directly.         c
c   2. If modifications are performed on the routines, these           c
c      modifications will be communicated to the authors.  The         c
c      modified routines will remain the sole property of the authors. c
c   3. Due acknowledgment must be made of the use of the FFSQP         c
c      routines in research reports or publications.  Whenever         c
c      such reports are released for public access, a copy should      c
c      be forwarded to the authors.                                    c
c   4. The FFSQP routines may only by used for research and            c
c      development, unless it has been agreed otherwise with the       c
c      authors in writing.                                             c
c                                                                      c
c Copyright (c) 1989 --- 1998 by Jian L. Zhou, Andre L. Tits,          c
c                                and Craig T. Lawrence                 c
c All Rights Reserved.                                                 c
c                                                                      c
c                                                                      c
c Enquiries should be directed to:                                     c
c                                                                      c
c      Prof. Andre L. Tits                                             c
c      Electrical Engineering Dept.                                    c
c      and Institute for Systems Research                              c
c      University of Maryland                                          c
c      College Park, Md 20742                                          c
c      U. S. A.                                                        c
c                                                                      c
c      Phone : 301-405-3669                                            c
c      Fax   : 301-405-6707                                            c
c      E-mail: andre@eng.umd.edu                                       c
c                                                                      c
c                                                                      c
c  References:                                                         c
c  [1] E. Panier and A. Tits, `On Combining Feasibility, Descent and   c
c      Superlinear Convergence In Inequality Constrained Optimization',c
c      Mathematical Programming 59(1993), 261-276.                     c
c  [2] J. F. Bonnans, E. Panier, A. Tits and J. Zhou, `Avoiding the    c
c      Maratos Effect by Means of a Nonmonotone Line search: II.       c
c      Inequality Problems - Feasible Iterates', SIAM J. Numer. Anal.  c
c      29(1992), 1187-1202.                                            c
c  [3] J.L. Zhou and A. Tits, `Nonmonotone Line Search for Minimax     c
c      Problems', J. Optim. Theory Appl.76(1993), 455-476.             c
c  [4] J.L. Zhou and A. Tits, `User's Guide for FFSQP Version 3.7 :    c
c      A Fortran Code for Solving Optimization Programs, Possibly      c
c      Minimax,with General Inequality Constraints and Linear Equality c
c      Constraints, Generating Feasible Iterates', Institute for       c
c      Systems Research, University of Maryland,Technical Report       c
c      SRC-TR-92-107r5, College Park, MD 20742, 1997.                  c
c  [5] C.T. Lawrence and A.L. Tits, `Nonlinear Equality Constraints    c
c      in Feasible Sequential Quadratic Programming,' Optimization     c
c      Methods and Software (6), March, 1996, pp. 265-282.             c
c                                                                      c
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
c
      integer i,io,ipp,iter,j,ncallg,ncallf,ncnstr,nclin,nctotl,leniw,
     *        lenw,nwx,nwbl,nwbu,nwgrg,nwgpf,nwpenp,nwa,nwcvec,nwhess,
     *        nwcla,nww,nrowa,modd,nppram,iwnob,iwncn,iwia,iwisp,iwist,
     *        iwiw,nwdi,nwd,nwff,nwgrf,nwclla,nwhes1,nwsp,nwbak,nwsg,M,
     *       maxit,nob,nobL,nnineq,info,idummy,nshift,max0,modem,lstype,
     *       nstop,initvl,nn,nnn,nwgm,ipsp,ipspan,ipyes,iprnto,mod
      double  precision epsmac,QLeps,small,xi,gi,gmax,dummy,big,tolfea,
     *        rteps,epskt,upert,valnom,dsqrt,dmax1
      logical feasbl,feasb,prnt,nspace,Linfty,nAD,rolis1,d0is0
      common  /fsqpp1/nnineq,M,ncallg,ncallf,modd,lstype,nstop,
     *        /fsqpp2/io,ipp,ipsp,ipyes,info,idummy,iter,initvl,
     *        /fsqpp3/epsmac,rteps,upert,valnom,
     *        /fsqpp4/rolis1,d0is0,
     *        /fsqpq1/big,tolfea,/fsqpq2/maxit
      common  /CMACHE/QLeps
c
c     compute the machine precision
c
      io=6
      rolis1=.false.
c     iwiw=6*nparam+8*max0(1,nineq+neq)+7*max0(nf,1)+30
c     i=nineq+neq+1
c     nww=4*nparam**2+5*i*nparam+3*(nf+1)*nparam+26*(nparam+nf+1)
c    *   +45*i+100
c     if(iwsize.ge.iwiw.and.nwsize.ge.nww) goto 10
c       if(iwsize.lt.iwiw) write(io,9906) iwiw
c       if(nwsize.lt.nww)  write(io,9907) nww
c       info=7
c       goto 9000
c
 10   iter=0
      nstop=1
      nn=nineqn+neqn
      epsmac=small()
      QLeps=epsmac
      tolfea=epsmac*1.d+02
      big=bigbnd
      rteps=dsqrt(epsmac)
      upert=udelta
c
      i=mod(iprint,10)
      ipspan=max0(iprint-i,1)
      iprnto=iprint
      if(iprint.ge.10) iprint=i
      if(iprint.lt.2) ipspan=1
      if(ipspan.lt.10) ipyes=0
      nob=0
      gmax=-bigbnd
      info=0
      ipsp=ipspan
      ipp=iprint
      ncnstr=nineq+neq
      nnineq=nineq
c
c     check input data
c
      if(iprint.gt.0) write(io,9900)
      call check(nparam,nf,Linfty,nAD,nineq,nineqn,neq,neqn,
     *           mode,modem,nwa,eps,bigbnd,bl,bu)
      if(info.eq.7) goto 9000
      lstype=nwa
c
      maxit=max0(max0(miter,10*max0(nparam,ncnstr)),1000)
      feasbl=.true.
      feasb=.true.
      prnt=.false.
      nspace=.false.
      nppram=nparam+1
      nshift=nparam**2+nppram**2
c
c  check whether x is within bounds
c
      do 100 i=1,nparam
        xi=x(i)
        if(bl(i).le.xi.and.bu(i).ge.xi) goto 100
        feasbl=.false.
        goto 110
 100  continue
 110  nclin=ncnstr-nn
c
c  check whether linear constraints are feasible
c
      if(nclin.eq.0) goto 210
      do 200 i=1,nclin
        j=i+nineqn
        if(j.le.nineq) then
          call constr(nparam,j,x,gi)
          if(gi.le.epsmac) goto120
          feasbl=.false.
        else if(j.gt.nineq) then
          call constr(nparam,j+neqn,x,gi)
          if(dabs(gi).le.epsmac) goto 120
          feasbl=.false.
        endif
 120    g(j)=gi
 200  continue
 210  if(feasbl) goto 240
      if(iprint.le.0) goto 230
        write(io,9901)
        call sbout1(io,nparam,'x                   ',dummy,x,2,1)
        prnt=.true.
 230  nctotl=nparam+nclin
      leniw=max0(2*nparam+2*nctotl+3,2*nclin+2*nparam+6)
      if(leniw.le.iwsize)then
        leniw=iwsize
      else
        write(io,9906) leniw
        info=7
        nspace=.true.
      endif
      nwx=1
      nwbl=nwx+nparam
      nwbu=nwbl+nctotl+4
      nwgrg=nwbu+nctotl+2
      nwa=nwgrg+nclin*nparam+1
      nwcvec=nwa+nparam*nclin+1
      nwhess=nwcvec+nparam
      nwcla=nwhess+nparam*nparam
      nww=nwcla+nctotl+nparam
      lenw=2*nparam**2+10*nparam+2*nctotl+1
      if((nww+lenw).le.nwsize) then
        lenw=nwsize-nww
        if(.not.nspace) goto 235
        write(io,9909)
        goto 9000
      else
        write (io,9907) nww+lenw
        write(io,9909)
        info=7
        goto 9000
      endif
c
c     attempt to generate a point satisfying all linear constraints
c
 235  nrowa=max0(nclin,1)
      call initpt(nparam,nineqn,neq,neqn,nclin,nctotl,nrowa,x,bl,bu,
     *            iw,leniw,w(nwx),w(nwbl),w(nwbu),g(nineqn+1),w(nwgrg),
     *            w(nwa),w(nwcvec),w(nwhess),w(nwcla),w(nwbl+nparam+3),
     *            w(nww),lenw,constr,gradcn)
      if(info.ne.0) goto 9000
 240  do 245 i=1, neq-neqn
 245    g(nineq+neqn+i)=g(nineq+i)
      if(nn.ne.0) goto 510
      goto 605
c
 290    do 300 i=1,nob
 300      w(i+nineqn+nshift)=w(i+nshift)
        nob=0
c
 510  continue
      if(info.eq.-1) goto 540
        do 520 i=1,nineqn
          call constr(nparam,i,x,w(i+nineqn+nshift))
          if(w(i+nineqn+nshift).gt.0.d0) feasb=.false.
 520    continue
        ncallg=nineqn
        if(feasb) goto 540
c
c     set indxob(i) in phase 1
c
        do 530 i=1,nineqn
          nob=nob+1
          iw(nob)=i
          w(nob+nshift)=w(i+nineqn+nshift)
          gmax=dmax1(gmax,w(nob+nshift))
 530    continue
      goto 580
 540    do 550 i=1,nineqn
          g(i)=w(i+nineqn+nshift)
          iw(nineqn+i+1)=i
 550    continue
        do 560 i=1,neq-neqn
          g(i+nineq+neqn)=g(i+nineq)
 560    continue
        do 570 i=1,neqn
          j=i+nineq
          call constr(nparam,j,x,g(j))
          iw(nineqn+nineqn+i+1)=j
 570    continue
        ncallg=ncallg+neqn
 580  continue
c
 605  if(iprint.le.0.or..not.feasb.or.prnt) goto 610
        write(io,9902)
        call sbout1(io,nparam,'x                   ',dummy,x,2,1)
        prnt=.true.
 610  if(nob.ne.0) goto 620
      if(iprint.le.0) goto 615
        if(info.eq.0) goto 613
        write(io,9904) ncallg
        if(ipp.eq.0) write(io,9910) iter
        if(ipp.gt.0) write(io,9910) iter-1
        if(ipp.eq.0) iter=iter+1
 613    if(.not.feasb.or.feasbl) goto 614
          write(io,9903)
          call sbout1(io,nparam,'x                   ',dummy,x,2,1)
 614    if(info.eq.0.and.prnt.and.feasb) goto 615
          write(io,9903)
          call sbout1(io,nparam,'x                   ',dummy,x,2,1)
 615  feasb=.true.
      feasbl=.true.
 620  nspace=.false.
      if(ipp.le.0.or.feasb.or.prnt) goto 630
        write(io,9901)
        call sbout1(io,nparam,'x                   ',dummy,x,2,1)
        prnt=.true.
 630  if(nob.eq.0) nob=1
c
c     set indxcn(1)--indxcn(ncnstr)
c
      if(feasb) nnn=nn
      if(.not.feasb) nnn=0
      do 700 i=1,nnn
 700    iw(nob+i)=iw(nineqn+i+1)
 710  do 800 i=1,nineq-nineqn
 800    iw(nob+nnn+i)=nineqn+i
      do 805 i=1,neq-neqn
        if(feasb) iw(nob+nineq+neqn+i)=nineq+neqn+i
        if(.not.feasb) iw(nineq+i)=nineq+neqn+i
 805  continue
      if(.not.feasb) goto 810
        nob=nf
        info=0
        ipp=iprint
        ipsp=ipspan
        modd=modem
        epskt=eps
        if(Linfty) nobL=2*nob
        if(.not.Linfty) nobL=nob
        if(nob.ne.0.or.neqn.ne.0) goto 910
        write(io,9908)
      goto 9000
 810    ipp=0
        ipsp=1
        modd=0
        nobL=nob
        info=-1
        epskt=1.d-10
 910  nctotl=nppram+ncnstr+max0(nobL,1)
      iwnob=1
      if(feasb) iwncn=iwnob+1
      if(.not.feasb) iwncn=iwnob+nob
      iwia=iwncn+ncnstr
      iwisp=iwia+nn+max0(nob,1)
      iwist=iwisp+nnineq-nineqn+1
      iwiw=iwist+nn+max0(nob,1)
      leniw=2*(ncnstr+max0(nobL,1))+2*nppram+6
c
      if((iwiw+leniw).le.iwsize) then
        leniw=iwsize-iwiw
      else
        write (io,9906) iwiw+leniw
        info=7
        nspace=.true.
      endif
      M=4
      if(modem.eq.1.and.nn.eq.0) M=3
      nwhess=1
      nwhes1=nwhess+nparam**2
      nwff=nwhes1+nppram**2
      nwx=nwff+max0(nob,1)+1
      nwdi=nwx+nppram
      nwd=nwdi+nppram
      nwgm=nwd+nppram
      nwgrg=nwgm+max0(1,4*neqn)
      nwgrf=nwgrg+ncnstr*nparam+1
      nwgpf=nwgrf+nparam*max0(nob,1)+1
      nwpenp=nwgpf+nparam
      nwa=nwpenp+neqn+1
      nwbl=nwa+(ncnstr+max0(nobL,1))*(nppram+1)
      nwbu=nwbl+nctotl+4
      nwcla=nwbu+nctotl+2
      nwclla=nwcla+nctotl+nppram
      nwcvec=nwclla+nctotl
      nwsp=nwcvec+nppram
      nwbak=nwsp+M+1
      nwsg=nwbak+max0(nob,1)+ncnstr+1
      nww=nwsg+neqn+1
      lenw=2*nppram*nppram+10*nppram+6*(ncnstr+max0(nobL,1)+1)
c
      if((nww+lenw).le.nwsize) then
        lenw=nwsize-nww
        if(.not.nspace) goto 920
        write(io,9909)
        goto 9000
      else
        write (io,9907) nww+lenw
        write(io,9909)
        info=7
        goto 9000
      endif
c
 920  do 1000 i=nwx,nwx+nparam-1
 1000   w(i)=x(i-nwx+1)
      w(nwx+nparam)=gmax
      if(.not.feasb) goto 1150
        do 1100 i=1,neqn
          if(g(i+nineq).gt.0d0) w(nwsg+i-1)=-1.d0
          if(g(i+nineq).le.0d0) w(nwsg+i-1)=1.d0
 1100   continue
c
c     either attempt to generate a point satisfying all constraints
c     or try to solve the original problem
c
 1150 nrowa=max0(ncnstr+max0(nobL,1),1)
      call FFSQP1(miter,nparam,nob,nobL,nineqn,neq,neqn,ncnstr,nctotl,
     *            nrowa,feasb,epskt,epseqn,bl,bu,iw(iwnob),iw(iwncn),
     *            iw(iwia),iw(iwisp),iw(iwist),iw(iwiw),leniw,w(nwx),
     *            w(nwdi),w(nwd),g,w(nwgm),w(nwgrg),w(nwff),w(nwgrf),
     *            w(nwgpf),w(nwpenp),w(nwa),w(nwbl),w(nwbu),w(nwcla),
     *            w(nwclla),w(nwcvec),w(nwbl+nparam+3),w(nwhess),
     *            w(nwhes1),w(nwsp),w(nwbak),w(nwsg),w(nww),lenw,
     *            obj,constr,gradob,gradcn)
      do 1200 i=1,nparam
 1200   x(i)=w(nwx+i-1)
      if(info.eq.-1) goto 290
      if(info.eq.0.or.feasb) goto 1220
        info=2
        write(io,9905)
      goto 9000
 1220 do 1300 i=1,nf
 1300   f(i)=w(nwff+i-1)
      if(nobL.le.1) idummy=0
      if(nobL.gt.1) idummy=1
      if(nf.eq.1) nob=0
      do 1400 i=1,nparam+ncnstr+nob
        j=i
        if(i.gt.nparam.and.i.le.(nparam+ncnstr))
     *    j=nparam+iw(iwncn+i-nparam)
        if(i.le.nparam) then
          w(i)=w(nwclla+j-1)
        else if(i.gt.nparam) then
          if(i.le.(nparam+ncnstr)) j=nparam+iw(iwncn+i-1-nparam)
          w(i)=w(nwclla+j-1+idummy)
        endif
 1400 continue
      if(nf.eq.1) w(nparam+ncnstr+1)=1.d0
c
 9000 inform=info
      iprint=iprnto
      return
 9900 format(1x,// 1x,'   FFSQP Version 3.7b  (Released January 1998)'
     *   /   1x,'            Copyright (c) 1989 --- 1998         '
     *   /   1x,'               J.L. Zhou, A.L. Tits,            '
     *   /   1x,'                 and C.T. Lawrence              '
     *   /   1x,'                All Rights Reserved             ',//)
 9901 format(1x,'The given initial point is infeasible for inequality',
     *       /10x,'constraints and linear equality constraints:')

 9902 format(1x,'The given initial point is feasible for inequality',
     *   /8x,'constraints and linear equality constraints:')
 9903 format(1x,'Starting from the generated point feasible for',
     *   ' inequality',
     *   /10x,'constraints and linear equality constraints:')
 9904 format(1x,'To generate a point feasible for nonlinear inequality',
     *   /1x,'constraints and linear equality constraints,',
     *   ' ncallg = ',i10)
 9905 format(1x,'Error: No feasible point is found for nonlinear',
     *   ' inequality',
     *    /8x,'constraints and linear equality constraints'/)
 9906 format(1x,'iwsize should be bigger than', i20)
 9907 format(1x,'nwsize should be bigger than', i20)
 9908 format(1x,'current feasible iterate with no objective specified'/)
 9909 format(1x,/)
 9910 format(43x,'iteration = ',i10)
      end
      subroutine FFSQP1(miter,nparam,nob,nobL,nineqn,neq,neqn,ncnstr,
     *                  nctotl,nrowa,feasb,epskt,epseqn,xl,xu,indxob,
     *                  indxcn,iact,iskip,istore,iw,leniw,x,di,d,g,gm,
     *                  gradg,f,gradf,grdpsf,penp,a,bl,bu,clamda,
     *                  cllamd,cvec,bj,hess,hess1,span,backup,signeq,
     *                  w,lenw,obj,constr,gradob,gradcn)
c
c     FFSQP : main routine for the optimization
c
c     implicit real*8(a-h,o-z)
      integer miter,nparam,nob,nobL,nineqn,neq,neqn,ncnstr,nctotl,nrowa,
     *        leniw,lenw
      integer indxob(1),indxcn(1),iact(1),iskip(1),
     *        istore(1),iw(leniw)
c     integer indxob(nob),indxcn(ncnstr),iact(nob+nineqn+neqn),iskip(1),
c    *        istore(nineqn+nob+neqn),iw(leniw)
      double  precision epskt,epseqn
      double  precision xl(nparam),xu(nparam),x(nparam+1),di(nparam+1),
     *        d(nparam+1),g(1),gm(1),gradg(nparam,1),
     *        f(1),gradf(nparam,1),grdpsf(nparam),penp(1),
     *        a(nrowa,1),bl(1),bu(1),clamda(1),
     *        cllamd(1),cvec(nparam+1),bj(1),
     *        hess(nparam,nparam),hess1(1),span(1),
     *        backup(1),signeq(1),w(lenw)
c     double  precision xl(nparam),xu(nparam),x(nparam+1),di(nparam+1),
c    *        d(nparam+1),g(ncnstr),gm(4*neqn),gradg(nparam,ncnstr),
c    *        f(nob),gradf(nparam,nob),grdpsf(nparam),penp(neqn),
c    *        a(nrowa,1),bl(nctotl),bu(nctotl),clamda(nctotl+nparam+1),
c    *        cllamd(nctotl),cvec(nparam+1),bj(nrowa),
c    *        hess(nparam,nparam),hess1(nparam+1,nparam+1),span(1),
c    *        backup(nob+ncnstr),signeq(neqn),w(lenw)
      external obj,constr,gradob,gradcn
      logical feasb
c
      integer nnineq,M,ncallg,ncallf,mode,io,iprint,info,ipd,iter,nstop,
     *        initvl,ipspan,ipyes,lstype
      double  precision bigbnd,tolfea,epsmac,rteps,udelta,valnom
      logical dlfeas,local,update,first,rolis1,d0is0
      common  /fsqpp1/nnineq,M,ncallg,ncallf,mode,lstype,nstop,
     *        /fsqpp2/io,iprint,ipspan,ipyes,info,ipd,iter,initvl,
     *        /fsqpp3/epsmac,rteps,udelta,valnom,
     *        /fsqpp4/rolis1,d0is0,
     *        /fsqpq1/bigbnd,tolfea,
c    *        /fsqp1/rentry,
     *        /fsqplo/dlfeas,local,update,first
c
c     bj(1+) is equivalent to bl(nparam+3+)
c
      integer i,ng,iskp,nfs,ncf,ncg,nn,non,nstart,nrst,ncnst1,nctot1
      double  precision Cbar,Ck,dbar,fmax,fM,fMp,steps,d0nm,dummy,
     *        sktnom,scvneq,grdftd,dmax1,psf
c
      initvl=1
      first=.true.
      nrst=0
      ipd=0
      if(nobL.gt.1) ng=1
      if(nobL.le.1) ng=0
      if(iter.eq.0) call diagnl(nparam,1.d0,hess)
      if(.not.feasb) goto 5
        first=.true.
        if(iter.gt.0) iter=iter-1
        if(iter.ne.0) call diagnl(nparam,1.d0,hess)
 5    Cbar=1.d-02
      Ck=Cbar
      dbar=5.0d0
      nstart=1
      ncallf=0
      nstop=1
      nfs=0
      non=miter
      if(mode.eq.0) goto 10
        nfs=M
        non=0
 10   if(feasb) then
        nn=nineqn+neqn
        ncnst1=ncnstr
        nctot1=nctotl
      else
        nn=0
        ncnst1=ncnstr-nineqn-neqn
        nctot1=nnineq-nineqn+neq-neqn+nparam
        if(nob.gt.1) nctot1=nctot1+1
      endif
      scvneq=0.d0
      do 100 i=1,ncnst1
        valnom=g(indxcn(i))
        backup(i)=valnom
        if(feasb.and.i.gt.nineqn.and.i.le.nn) then
          gm(i-nineqn)=valnom*signeq(i-nineqn)
          scvneq=scvneq+dabs(valnom)
        endif
        if(.not.feasb.or.i.gt.nn) goto 20
          iact(i)=indxcn(i)
          istore(i)=0
          if(i.gt.nineqn) penp(i-nineqn)=2.d0
 20     call gradcn(nparam,indxcn(i),x,gradg(1,indxcn(i)),constr)
 100  continue
      call nullvc(nparam,grdpsf)
      psf=0.d0
      if(.not.feasb.or.neqn.eq.0) goto 110
        call resign(nparam,neqn,psf,grdpsf,penp,g(nnineq+1),
     *              gradg(1,nnineq+1),signeq,12,12)
 110  fmax=-bigbnd
      if(feasb.and.nob.eq.0) then
        fmax=0.d0
        fMp=0.d0
      endif
      do 140 i=1,nob
        if(.not.feasb) goto 120
          iact(nn+i)=i
          istore(nn+i)=0
          call obj(nparam,i,x,f(i))
          valnom=f(i)
          backup(i+ncnst1)=valnom
          call gradob(nparam,i,x,gradf(1,i),obj)
          ncallf=ncallf+1
          if(nobL.ne.nob) fmax=dmax1(fmax,-f(i))
        goto 130
 120      valnom=f(i)
          iact(i)=i
          istore(i)=0
          call gradcn(nparam,indxob(i),x,gradf(1,i),constr)
 130    fmax=dmax1(fmax,f(i))
 140  continue
      fM=fmax
      fMp=fM-psf
      span(1)=fM
c
      if(iprint.lt.3.or..not.first.or.ipyes.gt.0) goto 600
        do 300 i=1,nob
          if(.not.feasb) goto 250
            if(nob.gt.1)
     *        call sbout2(io,nparam,i,'gradf(j,',')',gradf(1,i))
            if(nob.eq.1)
     *        call sbout1(io,nparam,'gradf(j)            ',
     *                    dummy,gradf(1,1),2,2)
          goto 300
 250        call sbout2(io,nparam,indxob(i),'gradg(j,',')',gradf(1,i))
 300    continue
 310    if(ncnstr.eq.0) goto 410
        do 400 i=1,ncnst1
 400      call sbout2(io,nparam,indxcn(i),'gradg(j,',')',
     *                gradg(1,indxcn(i)))
        if(neqn.eq.0) goto 410
        call sbout1(io,nparam,'grdpsf(j)           ',dummy,grdpsf,2,2)
        call sbout1(io,neqn,'P                   ',dummy,penp,2,2)
 410    do 500 i=1,nparam
 500      call sbout2(io,nparam,i,'hess (j,',')',hess(1,i))
c
c     main loop of the algorithm
c
 600  nstop=1
 601  continue
        call out(miter,nparam,nob,nineqn,nn,neqn,ncnst1,x,g,
     *           f,fmax,fM,psf,steps,sktnom,d0nm,feasb)
        if(nstop.ne.0) goto 810
        if(.not.feasb) goto 809
          do 700 i=1,ncnst1
 700        g(i)=backup(i)
          do 800 i=1,nob
 800        f(i)=backup(i+ncnst1)
          if(neqn.eq.0) goto 809
            do 805 i=1,neqn
              cllamd(nparam+ng+nnineq+i)=(cllamd(nparam+ng+nnineq+i)-
     *          penp(i))*signeq(i)
 805        continue
 809      return
 810    continue
        if(ipspan.ge.10.and.iprint.ge.2.and.ipyes.eq.0)
     *    write(io,9900) iter
        call dir(nparam,nob,nobL,nineqn,neq,neqn,nn,ncnst1,nctot1,nrowa,
     *           feasb,steps,epskt,epseqn,sktnom,scvneq,Ck,d0nm,grdftd,
     *           xl,xu,indxob,indxcn,iact,iskp,iskip,istore,iw,leniw,
     *           x,di,d,g,gradg,f,fmax,fM,fMp,psf,gradf,grdpsf,penp,a,
     *           bl,bu,clamda,cllamd,cvec,bj,hess,hess1,w,lenw,
     *           backup,signeq,obj,constr)
        if(nstop.eq.0) goto 601
        first=.false.
        if(update.or.d0is0) goto 820
        call step(nparam,nob,nobL,nineqn,neq,neqn,nn,ncnst1,ncg,ncf,
     *            indxob,indxcn,iact,iskp,iskip,istore,feasb,grdftd,
     *            f,fmax,fM,fMp,psf,penp,steps,scvneq,bu,x,di,d,g,w,
     *            backup,signeq,obj,constr)
        if(nstop.eq.0) goto 601
 820    call hesian(nparam,nob,nobL,nineqn,neq,neqn,nn,ncnst1,nctot1,
     *              nfs,nstart,feasb,bigbnd,bu,x,f,fmax,fM,fMp,psf,
     *              gradf,grdpsf,penp,g,gm,gradg,indxob,indxcn,cllamd,
     *              bl,clamda,di,hess,d,steps,nrst,signeq,span,
     *              obj,constr,gradob,gradcn,
     *              hess1,cvec,bj,w,lenw,iw,leniw)
        if(nstop.eq.0 .or. mode.eq.0) goto 601
        if(d0nm.gt.dbar) Ck=dmax1(dble(0.5*Ck),Cbar)
        if(d0nm.le.dbar.and.dlfeas) Ck=Ck
        if(d0nm.le.dbar.and..not.dlfeas.and..not.rolis1) Ck=10.0*Ck
      goto 601
 9900 format(1x,9hiteration,t22,i22)
      end
      subroutine check(nparam,nf,Linfty,nAD,nineq,nnl,neq,neqn,
     *                 mode,modem,lstype,eps,bigbnd,bl,bu)
c
c     FFSQP : check input data
c
c     implicit real*8(a-h,o-z)
      integer nparam,nf,nineq,nnl,neq,neqn,mode,modem,lstype
      double  precision bigbnd,eps
      double  precision bl(nparam),bu(nparam)
      logical Linfty,nAD
c
      integer io,iprint,ipspan,ipyes,info,idum1,idum2,idum3
      double  precision epsmac,dummy1,dummy2,dummy3
      common  /fsqpp2/io,iprint,ipspan,ipyes,info,idum1,idum2,idum3,
     *        /fsqpp3/epsmac,dummy1,dummy2,dummy3
c
      integer i
      double  precision bli,bui
c
      if (nparam.le.0)
     *  call error('nparam should be positive!              ',info,io)
      if (nf.lt.0)
     *  call error('nf     should not be negative!          ',info,io)
      if (nnl.lt.0)
     *  call error('nineqn should not be negative!          ',info,io)
      if (nineq.lt.nnl)
     *  call error('nineq  should be no smaller than nineqn!',info,io)
      if (neqn.lt.0)
     *  call error('neqn   should not be negative!          ',info,io)
      if (neq.lt.neqn)
     *  call error('neq    should not be smaller than neqn  ',info,io)
      if (nparam.le.(neq-neqn))
     *  call error('Need nparam >number of linear equalities',info,io)
      if (nparam.lt.neq) then
        call error('WARNING: nparam < neq                   ',info,io)
        info=0
      endif
      if (iprint.lt.0.or.iprint.gt.3)
     *  call error('iprint is not a valid number            ',info,io)
      if (eps.gt.epsmac) goto 10
      call error('eps    should be bigger than epsmac!    ',info,io)
      write(io,9902) epsmac
 10   if(mode.ge.200) then
        lstype=2
        mode=mode-100
      else
        lstype=1
      endif
      if (.not.(mode.eq.100.or.mode.eq.101.or.
     *          mode.eq.110.or.mode.eq.111))
     *  call error('mode   is not properly specified!       ',info,io)
      if (info.eq.0) goto 20
      write(io,9903)
      goto 9000
c
 20   do 30 i=1,nparam
        bli=bl(i)
        bui=bu(i)
        if (bli.le.bui) goto 25
        write(io,9901)
        info=7
 25     if (info.ne.0) goto 9000
        if (bli.lt.(-bigbnd)) bl(i)=-bigbnd
        if (bui.gt.bigbnd)    bu(i)=bigbnd
 30   continue
c
      i=mode-100
      if(i.lt.10) then
        modem=0
      else
        modem=1
        i=i-10
      endif
      if(i.eq.0) Linfty=.false.
      if(i.eq.1) Linfty=.true.
c
 9000 return
 9901 format(1x,'lower bounds should be smaller than upper bounds',/)
 9902 format(1x,'epsmac = ',e22.14,' which is machine dependent',/)
 9903 format(1x,'Error: Input parameters are not consistent',/)
      end
      subroutine initpt(nparam,nnl,neq,neqn,nclin,nctotl,nrowa,x0,
     *                  bndl,bndu,iw,leniw,x,bl,bu,g,gradg,a,cvec,hess,
     *                  clamda,bj,w,lenw,constr,gradcn)
c
c     FFSQP : generation of a feasible point satisfying
c             simple bounds and linear constraints
c
c     implicit real*8(a-h,o-z)
      integer nparam,nnl,neq,neqn,nclin,nctotl,nrowa,leniw,lenw
      integer iw(leniw)
      double  precision x0(nparam),bndl(nparam),bndu(nparam),x(nparam),
     *        bl(1),bu(1),g(1),gradg(nparam,1),
     *        a(nrowa,1),cvec(nparam),hess(nparam,nparam),
     *        clamda(1),bj(1),w(lenw)
c     double  precision x0(nparam),bndl(nparam),bndu(nparam),x(nparam),
c    *        bl(nctotl),bu(nctotl),g(nclin),gradg(nparam,nclin),
c    *        a(nrowa,nparam),cvec(nparam),hess(nparam,nparam),
c    *        clamda(nctotl+nparam),bj(nclin),w(lenw)
      external constr,gradcn
c
c     bj(1) is equivalent to bl(nparam+3)
c
      integer io,iprint,ipspan,ipyes,info,ipd,idum,idum2,maxit,
     *        nnineq,M,id2,id3,id4,id5,id6
      double  precision epsmac,rteps,udelta,valnom,big,tolfea,objeps,
     *        objrep,gLgeps
      logical xisnew
      common  /fsqpp1/nnineq,M,id2,id3,id4,id5,id6
     *        /fsqpp2/io,iprint,ipspan,ipyes,info,ipd,idum,idum2,
     *        /fsqpp3/epsmac,rteps,udelta,valnom,
     *        /fsqpq1/big,tolfea,/fsqpq2/maxit
      common  /fsqpus/objeps,objrep,gLgeps,xisnew
c
      integer i,j,infoql,mnn
      double  precision x0i
c
      info=1
      do 10 i=1,nclin
        valnom=g(i)
        j=i+nnl
        if(j.le.nnineq) call gradcn(nparam,j,x0,gradg(1,i),constr)
        if(j.gt.nnineq)
     *    call gradcn(nparam,j+neqn,x0,gradg(1,i),constr)
 10   continue
      do 20 i=1,nparam
        x0i=x0(i)
        bl(i)=bndl(i)-x0i
        bu(i)=bndu(i)-x0i
        cvec(i)=0.d0
 20   continue
      do 30 i=nclin,1,-1
 30     bj(nclin-i+1)=-g(i)
      do 60 i=nclin,1,-1
        do 50 j=1,nparam
 50       a(nclin-i+1,j)=-gradg(j,i)
 60   continue
      call diagnl(nparam,1.d0,hess)
      call nullvc(nparam,x)
C
      mnn=nrowa+2*nparam
      iw(1)=1
      call QL0001(nclin,neq-neqn,nrowa,nparam,nparam,mnn,hess,cvec,A,
     *            bj,bL,bU,X,clamda,io,infoql,0,w,lenw,iw,leniw)
      if(infoql.ne.0) goto 90
      do 70 i=1,nparam
 70     x0(i)=x0(i)+x(i)
      xisnew=.true.
      do 80 i=1,nclin
        j=i+nnl
        if(j.le.nnineq) call constr(nparam,j,x0,g(i))
        if(j.gt.nnineq) call constr(nparam,j+neqn,x0,g(i))
 80   continue
      info=0
 90   if(info.eq.1.and.iprint.ne.0) write(io,1000)
 1000 format(1x,'Error: No feasible point is found for the',
     *                 ' linear constraints',/)
      return
      end
      subroutine dir(nparam,nob,nobL,nineqn,neq,neqn,nn,ncnstr,nctotl,
     *               nrowa,feasb,steps,epskt,epseqn,sktnom,scvneq,Ck,
     *               d0nm,grdftd,xl,xu,indxob,indxcn,iact,iskp,iskip,
     *               istore,iw,leniw,x,di,d,g,gradg,f,fmax,fM,fMp,psf,
     *               gradf,grdpsf,penp,a,bl,bu,clamda,cllamd,cvec,bj,
     *               hess,hess1,w,lenw,backup,signeq,obj,constr)
c
c     FFSQP : computation of a search direction
c
c     implicit real*8(a-h,o-z)
      integer nparam,nob,nobL,nineqn,neq,neqn,nn,ncnstr,nctotl,nrowa,
     *        iskp,leniw,lenw
      integer indxob(1),indxcn(1),iact(1),iskip(1),
     *        istore(1),iw(leniw)
c     integer indxob(nob),indxcn(ncnstr),iact(nob+nineqn+neqn),iskip(1),
c    *        istore(nineqn+nob+neqn),iw(leniw)
      double  precision steps,epskt,epseqn,sktnom,Ck,d0nm,grdftd,
     *        fmax,fM,fMp,psf,scvneq
      double  precision xl(nparam),xu(nparam),x(nparam+1),di(nparam+1),
     *        d(nparam+1),g(1),gradg(nparam,1),f(1),
     *        gradf(nparam,1),grdpsf(nparam),penp(1),
     *        a(nrowa,nparam+1),bl(1),bu(1),clamda(1),cllamd(1),
     *        cvec(nparam+1),bj(nrowa),hess(nparam,nparam),
     *        hess1(nparam+1,nparam+1),w(lenw),
     *        backup(1),signeq(1)
c     double  precision xl(nparam),xu(nparam),x(nparam+1),di(nparam+1),
c    *        d(nparam+1),g(ncnstr),gradg(nparam,ncnstr),f(nob),
c    *        gradf(nparam,nob),
c    *        grdpsf(nparam),penp(neqn),a(nrowa,nparam+1),bl(nctotl),
c    *        bu(nctotl),clamda(nctotl+nparam+1),cllamd(nctotl),
c    *        cvec(nparam+1),bj(nrowa),hess(nparam,nparam),
c    *        hess1(nparam+1,nparam+1),w(lenw),
c    *        backup(nob+ncnstr),signeq(neqn)
      external obj,constr
      logical feasb
c
      integer nnineq,M,ncallg,ncallf,mode,io,iprint,ipspan,ipyes,info,
     *        ipd,iter,nstop,initvl,lstype
      double  precision epsmac,rteps,udelta,valnom,bigbnd,tolfea,
     *        objeps,objrep,gLgeps
      logical dlfeas,local,update,first,lqpsl,ld0,rolis1,d0is0,
     *        xisnew
      common  /fsqpp1/nnineq,M,ncallg,ncallf,mode,lstype,nstop,
     *        /fsqpp2/io,iprint,ipspan,ipyes,info,ipd,iter,initvl,
     *        /fsqpp3/epsmac,rteps,udelta,valnom
     *        /fsqpp4/rolis1,d0is0,
     *        /fsqpq1/bigbnd,tolfea,
     *        /fsqplo/dlfeas,local,update,first,
     *        /fsqpqp/lqpsl,ld0
      common  /fsqpus/objeps,objrep,gLgeps,xisnew
c
c     bj(1) is equivalent to bl(nparam+3)
c
      integer i,j,k,kk,ncg,ncf,nqprm0,nclin0,nctot0,infoqp,nqprm1,ncl,
     *        nclin1,nctot1,ncc,nff,nrowa0,nrowa1,ninq,nobb,nobbL,nncn,
     *        max0
      double  precision fmxl,vv,dx,dmx,dnm1,dnm,v0,v1,vk,temp1,temp2,
     *        theta,rhol,rhog,rho,grdfd0,grdfd1,dummy,grdgd0,grdgd1,
     *        thrshd,sign,scaprd,slope,lfuscp,dsqrt,dmin1,dmax1,dabs,
     *        adummy(1),dnmtil
      logical ltem1,ltem2,needd1
c
      ncg=0
      ncf=0
      iskp=0
      ncl=nnineq-nineqn
      local=.false.
      update=.false.
      lqpsl=.false.
      thrshd=tolfea
      needd1=.true.
      rolis1=.false.
c
      if(nobL.le.1) goto 10
        nqprm0=nparam+1
        nclin0=ncnstr+nobL
      goto 20
 10     nqprm0=nparam
        nclin0=ncnstr
 20   nctot0=nqprm0+nclin0
      vv=0.d0
      nrowa0=max0(nclin0,1)
      do 25 i=1,ncnstr
        if(feasb) then
          if(i.gt.nineqn.and.i.le.nnineq) iskip(nnineq+2-i)=i
          iw(i)=i
        else if(.not.feasb) then
          if(i.le.ncl) iskip(ncl+2-i)=nineqn+i
          if(i.le.ncl) iw(i)=nineqn+i
          if(i.gt.ncl) iw(i)=nineqn+neqn+i
        endif
 25   continue
      do 27 i=1,nob
 27     iw(ncnstr+i)=i
      ld0=.true.
      call nullvc(nparam,cvec)
      d0is0=.false.
      call dqp(nparam,nqprm0,nob,nobL,nineqn,neq,neqn,nn,ncnstr,nclin0,
     *         nctot0,nrowa0,infoqp,iw,leniw,x,di,xl,xu,feasb,f,fmax,
     *         gradf,grdpsf,g,gradg,a,cvec,bl,bu,clamda,cllamd,bj,
     *         hess,hess1,di,w,lenw,vv,0)
      ld0=.false.
      if(infoqp.eq.0) goto 30
        info=5
        if(.not.feasb) info=2
        nstop=0
        goto 9000
c
c    reorder indexes of constraints and objectives
c
 30   if(nn.le.1) goto 45
      j=1
      k=nn
      do 40 i=nn,1,-1
        if(lfuscp(cllamd(nqprm0+indxcn(i)),thrshd).ne.0) then
          iact(j)=indxcn(i)
          j=j+1
        else
          iact(k)=indxcn(i)
          k=k-1
        endif
 40   continue
 45   if(nobL.le.1) goto 60
      j=nn+1
      k=nn+nob
      do 50 i=nob,1,-1
        kk=nqprm0+ncnstr
        ltem1=lfuscp(cllamd(kk+i),thrshd).ne.0
        ltem2=nobL.ne.nob.and.(lfuscp(cllamd(kk+i+nob),thrshd).ne.0)
        if(ltem1.or.ltem2) then
          iact(j)=i
          j=j+1
        else
          iact(k)=i
          k=k-1
        endif
 50   continue
c
 60   vv=f(iact(nn+1))
      d0nm=dsqrt(scaprd(nparam,di,di))
      if(.not.first.or.nclin0.ne.0) goto 110
        dx=dsqrt(scaprd(nparam,x,x))
        dmx=dmax1(dx,1.d0)
        if(d0nm.le.dmx) goto 110
        do 100 i=1,nparam
 100      di(i)=di(i)*dmx/d0nm
        d0nm=dmx
 110  call matrvc(nparam,nparam,nparam,nparam,hess,di,w)
      if(nn.eq.0) grdftd=-scaprd(nparam,w,di)
      sktnom=dsqrt(scaprd(nparam,w,w))
      if(gLgeps.gt.0.d0.and.sktnom.le.gLgeps) goto 115
      if(d0nm.gt.epskt) goto 120
 115  if(neqn.ne.0.and.scvneq.gt.epseqn) goto 120
        nstop=0
        if(.not.feasb) info=2
        if(iprint.lt.3.or.ipyes.gt.0) goto 9000
        if(nobL.le.1) nff=1
        if(nobL.gt.1) nff=2
        call sbout1(io,nparam,'multipliers  for  x ',dummy,cllamd,2,2)
        if(ncnstr.ne.0) call sbout1(io,ncnstr,'             for  g ',
     *                              dummy,cllamd(nparam+nff),2,2)
        if(nobL.gt.1) call sbout1(io,nob,'             for  f ',
     *                            dummy,cllamd(nparam+nff+ncnstr),2,2)
        goto 9000
 120  if(iprint.lt.3.or.ipyes.gt.0) goto 125
        call sbout1(io,nparam,'d0                  ',dummy,di,2,2)
        call sbout1(io,0,'d0norm              ',d0nm,adummy,1,2)
        call sbout1(io,0,'ktnorm              ',sktnom,adummy,1,2)
 125  temp1=dmin1(0.5d0*epskt,0.1d-1*rteps)
      if(neqn.eq.0.or.scvneq.le.epseqn.or.d0nm.gt.temp1) goto 127
        d0is0=.true.
        goto 9000
c
c     single objective without nonlinear constraints requires
c     no d1 and dtilde; multi-objectives without nonlinear
c     constraints requires no d1
c
 127  call nullvc(nparam,w)
      if(nn.ne.0) grdftd=slope(nob,nobL,neqn,nparam,feasb,f,gradf,
     *                         grdpsf,di,w,fmax,dummy,0)
      if(nn.eq.0.and.nobL.le.1) goto 1130
      if(nn.ne.0) goto 130
        dnm=d0nm
        rho=0.d0
        rhog=0.d0
        goto 310
c
c     compute modified first order direction d1
c
 130  if (mode.eq.1) then
        needd1=.false.
        vk=dmin1(Ck*d0nm*d0nm,d0nm)
        do 131 i=1,nn
          grdgd0=scaprd(nparam,gradg(1,indxcn(i)),di)
          temp1=vk+g(indxcn(i))+grdgd0
          if(temp1.le.0.d0) goto 131
          needd1=.true.
          goto 132
 131    continue
      endif
 132  if (needd1) then
        nqprm1=nparam+1
        if(mode.eq.0) nclin1=ncnstr+max0(nobL,1)
        if(mode.eq.1) nclin1=ncnstr
        nctot1=nqprm1+nclin1
        nrowa1=max0(nclin1,1)
        ninq=nnineq
        call di1(nparam,nqprm1,nob,nobL,nineqn,neq,neqn,ncnstr,nclin1,
     *           nctot1,nrowa1,infoqp,mode,iw,leniw,x,di,xl,xu,f,fmax,
     *           gradf,grdpsf,g,gradg,cvec,a,bl,bu,clamda,bj,hess1,d,
     *           w,lenw)
        if(infoqp.eq.0) goto 140
          info=6
          if(.not.feasb) info=2
          nstop=0
          goto 9000
 140    dnm1=dsqrt(scaprd(nparam,d,d))
        if(iprint.lt.3.or.ipyes.gt.0) goto 145
          call sbout1(io,nparam,'d1                  ',dummy,d,2,2)
          call sbout1(io,0,'d1norm              ',dnm1,adummy,1,2)
      else
        dnm1=0.d0
        do 141 i=1,nparam
 141      d(i)=0.d0
      endif
 145  if(mode.eq.1) goto 150
        v0=d0nm**2.1
        v1=dmax1(dble(0.5),dble(dnm1**2.5))
        rho=v0/(v0+v1)
        rhog=rho
      goto 250
 150    rhol=0.d0
        if(.not.needd1) goto 210
        do 200 i=1,nn
          grdgd0=scaprd(nparam,gradg(1,indxcn(i)),di)
          grdgd1=scaprd(nparam,gradg(1,indxcn(i)),d)
          temp1=vk+g(indxcn(i))+grdgd0
          temp2=grdgd1-grdgd0
          if(temp1.le.0.d0) goto 200
          if(temp2.ge.0.d0) goto 190
          rhol=dmax1(rhol,-temp1/temp2)
          if(rhol.lt.1.d0) goto 200
 190        rhol=1.0d0
            rolis1=.true.
            goto 210
 200    continue
 210    theta=0.2d0
        if(rhol.ne.0.d0) goto 220
c
c       to check if rhol is reset
c
          rhog=0.d0
          rho=0.d0
          dnm=d0nm
        goto 310
 220    if(nobL.gt.1) goto 230
          grdfd0=grdftd
          if(nob.eq.1) grdfd1=scaprd(nparam,gradf(1,1),d)
          grdfd1=grdfd1-scaprd(nparam,grdpsf,d)
          temp1=grdfd1-grdfd0
          if(temp1.le.0.d0) then
            rhog=rhol
          else
            rhog=dmin1(rhol,(theta-1.d0)*grdfd0/temp1)
          endif
        goto 240
 230      rhog=slope(nob,nobL,neqn,nparam,feasb,f,gradf(1,1),grdpsf,
     *               di,d,fmax,theta,mode)
          rhog=dmin1(rhol,rhog)
 240    rho=rhog
        if (steps.eq.1.d0.and.rhol.lt.0.5d0) rho=rhol
 250  continue
      do 300 i=1,nparam
        if (rho.ne.rhog) cvec(i)=di(i)
        di(i)=(1.d0-rho)*di(i)+rho*d(i)
 300  continue
      dnm=dsqrt(scaprd(nparam,di,di))
      if(iprint.lt.3.or.mode.eq.1.or.nn.eq.0.or.ipyes.gt.0) goto 310
        call sbout1(io,0,'rho                 ',rho,adummy,1,2)
        call sbout1(io,nparam,'d                   ',dummy,di,2,2)
        call sbout1(io,0,'dnorm               ',dnm,adummy,1,2)
 310  continue
 320  do 400 i=1,nob
 400    bl(i)=f(i)
      if (rho.eq.1.d0) goto 510
      if(nn.eq.0.or.iprint.ne.3.or.mode.eq.0.or.ipyes.gt.0) goto 410
        call sbout1(io,0,'Ck                  ',Ck,adummy,1,2)
        call sbout1(io,0,'rhol                ',rho,adummy,1,2)
        call sbout1(io,nparam,'dl                  ',dummy,di,2,2)
        call sbout1(io,0,'dlnorm              ',dnm,adummy,1,2)
 410  if(mode.eq.0) goto 510
        local=.true.
        call step(nparam,nob,nobL,nineqn,neq,neqn,nn,ncnstr,ncg,ncf,
     *            indxob,indxcn,iact,iskp,iskip,istore,feasb,grdftd,
     *            f,fmax,fM,fMp,psf,penp,steps,scvneq,bu,x,di,d,g,w,
     *            backup,signeq,obj,constr)
        if(update) goto 9000
        local=.false.
        if(rho.eq.rhog.or.nn.eq.0) goto 510
        do 500 i=1,nparam
 500      di(i)=(1-rhog)*cvec(i)+rhog*d(i)
        dnm=dsqrt(scaprd(nparam,di,di))
 510  if (nn.eq.0.or.iprint.lt.3.or.mode.eq.0.or.ipyes.gt.0) goto 520
        call sbout1(io,0,'rhog                ',rhog,adummy,1,2)
        call sbout1(io,nparam,'dg                  ',dummy,di,2,2)
        call sbout1(io,0,'dgnorm              ',dnm,adummy,1,2)
 520  if(rho.ne.0.d0) grdftd=slope(nob,nobL,neqn,nparam,feasb,bl,
     *                             gradf,grdpsf,di,d,fmax,theta,0)
      if(mode.eq.1.and.rho.eq.rhog) goto 610
      do 600 i=1,nparam
 600    bu(i)=x(i)+di(i)
      xisnew=.true.
 610  if(rho.ne.rhog) ncg=0
      ncc=ncg+1
      fmxl=-bigbnd
      ninq=ncg
      nncn=ncg
      j=0
c
c     iskip(1) --- iskip(iskp) store the indexes of linear inequality
c     constraints that are not to be used to compute d~
c     iskip(nnineq-nineqn+1) --- iskip(nnineq-ncn+1-iskp) store those
c     that are to be used to compute d~
c
      do 700 i=ncc,ncnstr
        if(i.le.nn) then
          kk=iact(i)
        else
          kk=indxcn(i)
        endif
        if(kk.le.nineqn.or.kk.gt.nnineq) goto 615
          iskip(ncl+1-j)=kk
          j=j+1
 615    if(kk.gt.nnineq) goto 617
        temp1=-0.2d0*(dnm*dsqrt(scaprd(nparam,gradg(1,kk),gradg(1,kk))))
        temp2=cllamd(nqprm0+kk)
        if(temp2.eq.0.d0.and.g(kk).lt.temp1) goto 620
 617      ninq=ninq+1
          iw(ninq)=kk
          if(feasb.and.kk.le.nineqn) istore(kk)=1
          call constr(nparam,kk,bu,g(kk))
          if(.not.feasb.or.feasb.and.kk.gt.(nnineq+neqn)) goto 700
          if(kk.le.nineqn) nncn=ninq
          fmxl=dmax1(fmxl,g(kk))
          if(.not.feasb) goto 618
          if(kk.le.nineqn.or.kk.gt.nnineq.and.kk.le.(nnineq+neqn))
     *       ncallg=ncallg+1
 618      if(dabs(fmxl).gt.bigbnd) goto 1130
        goto 700
 620      if(kk.le.nineqn) goto 700
          iskp=iskp+1
          iskip(iskp)=kk
          j=j-1
 700  continue
      if(neqn.ne.0) call resign(nparam,neqn,psf,grdpsf,penp,g(nnineq+1),
     *                          gradg(1,nnineq+1),signeq,10,20)
      ninq=ninq-neq
      if(.not.feasb) ninq=ninq+neqn
      if(ncg.eq.0) goto 810
      do 800 i=1,ncg
        iw(i)=iact(i)
        if(iact(i).le.nineqn) istore(iact(i))=1
        fmxl=dmax1(fmxl,g(iact(i)))
        if(dabs(fmxl).gt.bigbnd) goto 1130
 800  continue
 810  if(nobL.gt.1) goto 820
        iw(1+ninq+neq)=1
        nobb=nob
        goto 1110
 820  if(rho.ne.rhog) ncf=0
      nff=ncf+1
      nobb=ncf
      sign=1.d0
      fmxl=-bigbnd
      if(cllamd(nqprm0+ncnstr+iact(nn+1)).lt.0.d0) sign=-1.d0
      do 1000 i=nff,nob
        kk=iact(nn+i)
        if(.not.feasb) kk=iact(i)
        if(feasb) k=nn+1
        if(.not.feasb) k=1
        do 900 j=1,nparam
 900      w(j)=sign*gradf(j,iact(k))-gradf(j,kk)
        temp1=dabs(f(kk)-sign*vv)
        temp2=dnm*dsqrt(scaprd(nparam,w,w))
        if(temp1.eq.0.d0.or.temp2.eq.0.d0) goto 910
        temp1=temp1/temp2
        temp2=cllamd(nqprm0+ncnstr+kk)
        if(temp2.eq.0.d0.and.temp1.gt.0.2d0) goto 1000
 910    nobb=nobb+1
        iw(nobb+ninq+neq)=kk
        if(feasb)      istore(nineqn+kk)=1
        if(.not.feasb) istore(kk)=1
        if(.not.feasb) goto 920
          call obj(nparam,kk,bu,f(kk))
          ncallf=ncallf+1
          if(nobL.ne.nob) fmxl=dmax1(fmxl,-f(kk))
        goto 930
 920      call constr(nparam,indxob(kk),bu,f(kk))
          ncallg=ncallg+1
 930    fmxl=dmax1(fmxl,f(kk))
        if(dabs(fmxl).gt.bigbnd) goto 1130
 1000 continue
      if(ncf.eq.0) goto 1110
      do 1100 i=1,ncf
        iw(ninq+neq+i)=iact(i+nn)
        istore(nineqn+iact(i+nn))=1
        fmxl=dmax1(fmxl,f(iact(i+nn)))
        if(nobL.ne.nob) fmxl=dmax1(fmxl,-f(iact(i+nn)))
        if(dabs(fmxl).gt.bigbnd) goto 1130
 1100 continue
 1110 call matrvc(nparam,nparam,nparam,nparam,hess,di,cvec)
      vv=-dmin1(0.01d0*dnm,dnm**2.5)
c
c     compute a correction dtilde to d=(1-rho)d0+rho*d1
c
      if(nobL.ne.nob) nobbL=2*nobb
      if(nobL.eq.nob) nobbL=nobb
      if(nobbL.le.1) goto 1115
        nqprm0=nparam+1
        nclin0=ninq+neq+nobbL
      goto 1117
 1115   nqprm0=nparam
        nclin0=ninq+neq
 1117 nctot0=nqprm0+nclin0
      nrowa0=max0(nclin0,1)
      i=ninq+neq
      call dqp(nparam,nqprm0,nobb,nobbL,nncn,neq,neqn,nn,i,nclin0,
     *         nctot0,nrowa0,infoqp,iw,leniw,x,di,xl,xu,feasb,f,fmxl,
     *         gradf,grdpsf,g,gradg,a,cvec,bl,bu,clamda,cllamd,bj,
     *         hess,hess1,d,w,lenw,vv,1)
      if(infoqp.ne.0) goto 1130
      dnmtil=dsqrt(scaprd(nparam,d,d))
      if(dnmtil.gt.dnm) goto 1130
      if(dnmtil.eq.0.d0) goto 1119
        do 1118 i=1,nineqn+nob
 1118     istore(i)=0
 1119 if(iprint.lt.3.or.ipyes.gt.0) goto 9000
        call sbout1(io,nparam,'dtilde              ',dummy,d,2,2)
        call sbout1(io,0,'dtnorm              ',dnmtil,adummy,1,2)
        goto 9000
c
 1130 do 1200 i=1,nparam
 1200   d(i)=0.d0
      dnmtil=0.d0
 9000 return
      end
c
      subroutine dqp(nparam,nqpram,nob,nobL,nineqn,neq,neqn,nn,ncnstr,
     *               nclin,nctotl,nrowa,infoqp,iw,leniw,x0,di,xl,xu,
     *               feasb,f,fmax,gradf,grdpsf,g,gradg,a,cvec,bl,bu,
     *               clamda,cllamd,bj,hess,hess1,x,w,lenw,vv,job)
c     implicit double precision(a-h,o-z)
      integer nparam,nqpram,nob,nobL,nineqn,neq,neqn,nn,ncnstr,nclin,
     *        nctotl,nrowa,infoqp,leniw,lenw,job
      integer iw(leniw)
      double  precision fmax,vv
      double  precision x0(nparam),di(1),xl(nparam),xu(nparam),
     *        f(1),gradf(nparam,1),grdpsf(nparam),g(1),
     *        gradg(nparam,1),
     *        a(nrowa,1),cvec(1),bl(1),bu(1),clamda(1),
     *        cllamd(1),bj(1),hess(nparam,nparam),
     *        hess1(nparam+1,nparam+1),x(1),w(lenw)
c     double  precision x0(nparam),di(nqpram),xl(nparam),xu(nparam),
c    *        f(nob),gradf(nparam,nob),grdpsf(nparam),g(ncnstr),
c    *        gradg(nparam,ncnstr),
c    *        a(nrowa,nqpram),cvec(nqpram),bl(nctotl),bu(nctotl),
c    *        clamda(nctotl+nqpram),cllamd(nctotl),bj(nrowa),
c    *        hess(nparam,nparam),hess1(nparam+1,nparam+1),
c    *        x(nqpram),w(lenw)
      logical feasb
c
      integer io,idum1,idum2,idum3,idum4,idum5,idum6,idum7
      double  precision bigbnd,dummy,epsmac,rteps,dummy1,dummy2
      common  /fsqpp2/io,idum1,idum2,idum3,idum4,idum5,idum6,idum7,
     *        /fsqpp3/epsmac,rteps,dummy1,dummy2,
     *        /fsqpq1/bigbnd,dummy
c
c     bj(1) is equivalent to bl(nparam+3)
c
c     job=0 : compute d0; job=1 : compute  d~
c
      integer i,ii,j,iout,mnn,nqnp
      double  precision x0i,xdi
c
      iout=io
      do 100 i=1,nparam
        x0i=x0(i)
        if(job.eq.1) xdi=di(i)
        if(job.eq.0) xdi=0.d0
        bl(i)=xl(i)-x0i-xdi
        bu(i)=xu(i)-x0i-xdi
        cvec(i)=cvec(i)-grdpsf(i)
 100  continue
      if(nobL.le.1) goto 110
        bl(nqpram)=-bigbnd
        bu(nqpram)=bigbnd
 110  ii=ncnstr-nn
c
c     constraints are assigned to a in reverse order
c
      do 300 i=1,ncnstr
        x0i=vv
        if(i.le.(neq-neqn).or.(i.gt.neq.and.i.le.(ncnstr-nineqn)))
     *    x0i=0.d0
        if(.not.feasb) x0i=0.d0
        bj(i)=x0i-g(iw(ncnstr+1-i))
        do 200 j=1,nparam
 200      a(i,j)=-gradg(j,iw(ncnstr+1-i))
        if(nobL.gt.1) a(i,nqpram)=0.d0
 300  continue
      if(nobL.le.1) goto 510
      do 500 i=1,nob
        ii=ncnstr+i
        bj(ii)=fmax-f(iw(ii))
        if(nobL.gt.nob) bj(ii+nob)=fmax+f(iw(ii))
        do 400 j=1,nparam
          a(ii,j)=-gradf(j,iw(ii))
          if(nobL.gt.nob) a(ii+nob,j)=gradf(j,iw(ii))
 400    continue
        a(ii,nqpram)=1.d0
        if(nobL.gt.nob) a(ii+nob,nqpram)=1.d0
 500  continue
      cvec(nqpram)=1.d0
      goto 610
 510  if(nob.eq.0) goto 610
      do 600 i=1,nparam
 600    cvec(i)=cvec(i)+gradf(i,1)
 610  call matrcp(nparam,hess,nparam+1,hess1)
      call nullvc(nqpram,x)
c
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
c  The following modification is done inside QP0001
c  for the ease of interfacing with QPSOL
c
c     if(hess1(nqpram,nqpram).lt.qleps) hess1(nqpram,nqpram)=qleps
C
      iw(1)=1
      mnn=nclin+2*nqpram
      call QL0001(nclin,neq-neqn,nrowa,nqpram,nparam+1,mnn,hess1,cvec,A,
     *            bj,bL,bU,X,clamda,iout,infoqp,0,w,lenw,iw,leniw)
C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
      if(infoqp.ne.0.or.job.eq.1) goto 9000
      do 700 i=1,nqpram
        ii=nclin+i
        if(clamda(ii).eq.0.d0.and.clamda(ii+nqpram).eq.0.d0) then
          goto 700
        else if(clamda(ii).ne.0.d0) then
          clamda(ii)=-clamda(ii)
        else
          clamda(ii)=clamda(ii+nqpram)
        endif
 700  continue
      nqnp=nqpram+ncnstr
      do 800 i=1,nctotl
        if(i.le.nqpram) then
          ii=nclin+i
        else if(i.gt.nqpram.and.i.le.nqnp) then
          ii=nqnp+1-i
        else if(i.gt.nqnp) then
          ii=i-nqpram
        endif
        cllamd(i)=clamda(ii)
 800  continue
      if(nobL.eq.nob) goto 9000
      do 900 i=1,nob
        ii=i+nqpram+ncnstr
        cllamd(ii)=cllamd(ii)-cllamd(ii+nob)
 900  continue
 9000 return
      end
c
      subroutine di1(nparam,nqpram,nob,nobL,nineqn,neq,neqn,ncnstr,
     *               nclin,nctotl,nrowa,infoqp,mode,iw,leniw,x0,d0,
     *               xl,xu,f,fmax,gradf,grdpsf,g,gradg,cvec,a,bl,bu,
     *               clamda,bj,hess1,x,w,lenw)
c     implicit real*8(a-h,o-z)
      integer nparam,nqpram,nob,nobL,nineqn,neq,neqn,ncnstr,nclin,
     *        nctotl,nrowa,infoqp,mode,leniw,lenw,iw(leniw)
      double  precision fmax
      double  precision x0(nparam),d0(nparam),xl(nparam),xu(nparam),
     *        f(1),gradf(nparam,1),grdpsf(nparam),g(1),
     *        gradg(nparam,1),cvec(1),a(nrowa,1),
     *        bl(1),bu(1),clamda(1),bj(1),
     *        hess1(nparam+1,nparam+1),x(1),w(lenw)
c     double  precision x0(nparam),d0(nparam),xl(nparam),xu(nparam),
c    *        f(nob),gradf(nparam,nob+1),grdpsf(nparam),g(ncnstr),
c    *        gradg(nparam,ncnstr),cvec(nqpram),a(nrowa,nqpram),
c    *        bl(nctotl),bu(nctotl),clamda(nctotl+nqpram),bj(nrowa),
c    *        hess1(nparam+1,nparam+1),x(nqpram),w(lenw)
c
      integer io,idum1,idum2,idum3,idum4,idum5,idum6,idum7
      double  precision epsmac,rteps,dumm1,dumm2,bigbnd,dummy
      common  /fsqpp2/io,idum1,idum2,idum3,idum4,idum5,idum6,idum7,
     *        /fsqpp3/epsmac,rteps,dumm1,dumm2,
     *        /fsqpq1/bigbnd,dummy
c
c     bj(1) is equivalent to bl(nparam+3)
c
      integer i,ii,iout,j,mnn
      double  precision x0i,eta
c
      iout=io
      if(mode.eq.0) eta=0.1d0
      if(mode.eq.1) eta=3.d0
      do 100 i=1,nparam
        x0i=x0(i)
        bl(i)=xl(i)-x0i
        bu(i)=xu(i)-x0i
        if(mode.eq.0) cvec(i)=-eta*d0(i)
        if(mode.eq.1) cvec(i)=0.d0
 100  continue
      bl(nqpram)=-bigbnd
      bu(nqpram)=bigbnd
      cvec(nqpram)=1.d0
      ii=ncnstr-nineqn
      do 400 i=1,ncnstr
        bj(i)=-g(ncnstr+1-i)
        do 300 j=1,nparam
 300      a(i,j)=-gradg(j,ncnstr+1-i)
        a(i,nqpram)=0.d0
        if((i.gt.(neq-neqn).and.i.le.neq).or.i.gt.ii) a(i,nqpram)=1.d0
 400  continue
      if(mode.eq.1) goto 610
      i=0
 450  continue
        i=i+1
        ii=ncnstr+i
        if(nob.eq.0) bj(ii)=fmax
        if(nob.gt.0) bj(ii)=fmax-f(i)
        do 500 j=1,nparam
          if(nob.eq.0) a(ii,j)=grdpsf(j)
          if(nob.gt.0) a(ii,j)=-gradf(j,i)+grdpsf(j)
          if(nobL.gt.nob) a(ii+nob,j)=gradf(j,i)+grdpsf(j)
 500    continue
        a(ii,nqpram)=1.d0
        if(nobL.gt.nob) a(ii+nob,nqpram)=1.d0
        if(i.lt.nob) goto 450
 610  call diagnl(nqpram,eta,hess1)
      call nullvc(nqpram,x)
      hess1(nqpram,nqpram)=0.d0
c
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  The following modification is done inside QP0001
c  for the ease of interfacing with QPSOL
c
c     hess1(nqpram,nqpram)=qleps
C
      mnn=nclin+2*nqpram
      iw(1)=1
      call QL0001(nclin,neq-neqn,nrowa,nqpram,nparam+1,mnn,hess1,cvec,A,
     *             bj,bL,bU,X,clamda,iout,infoqp,0,w,lenw,iw,leniw)
C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      return
      end

      subroutine step(nparam,nob,nobL,nineqn,neq,neqn,nn,ncnstr,ncg,
     *                ncf,indxob,indxcn,iact,iskp,iskip,istore,feasb,
     *                grdftd,f,fmax,fM,fMp,psf,penp,steps,scvneq,xnew,
     *                x,di,d,g,w,backup,signeq,obj,constr)
c
c     FFSQP : Armijo or nonmonotone line search, with
c             some ad hoc strategies to decrease the number
c             of function evaluation as much as possible
c
c     implicit real*8(a-h,o-z)
      integer nparam,nob,nobL,nineqn,neq,neqn,nn,ncnstr,ncg,ncf,iskp
      integer indxob(1),indxcn(1),iact(1),iskip(1),
     *        istore(1)
c     integer indxob(nob),indxcn(ncnstr),iact(nn+nob),iskip(4),
c    *        istore(nineqn+nob)
      double  precision grdftd,fmax,fM,fMp,steps,scvneq,psf
      double  precision xnew(nparam),x(nparam),di(nparam),d(nparam),
     *        f(1),penp(1),g(1),w(1),backup(1),
     *        signeq(1)
c     double  precision xnew(nparam),x(nparam),di(nparam),d(nparam),
c    *        f(nob),penp(neqn),g(ncnstr),w(1),backup(nob+ncnstr),
c    *        signeq(neqn)
      external obj,constr
      logical feasb
c
      integer nnineq,M,ncallg,ncallf,mode,io,iprint,ipspan,ipyes,info,
     *        idum1,idum2,idum3,nstop,lstype
      double  precision epsmac,bigbnd,tolfea,dum1,dum2,dum3,fii,
     *        objeps,objrep,gLgeps
      logical lqpsl,ldummy,dlfeas,local,update,ldumm2,xisnew
      common  /fsqpp1/nnineq,M,ncallg,ncallf,mode,lstype,nstop,
     *        /fsqpp2/io,iprint,ipspan,ipyes,info,idum1,idum2,idum3,
     *        /fsqpp3/epsmac,dum1,dum2,dum3,
     *        /fsqpq1/bigbnd,tolfea,
     *        /fsqplo/dlfeas,local,update,ldumm2,
     *        /fsqpqp/lqpsl,ldummy
      common  /fsqpus/objeps,objrep,gLgeps,xisnew
c
      integer i,ii,ij,itry,ikeep,j,job,nlin,mnm,ntwo
      double  precision prod1,prod,dummy,tolfe,dmax1,ostep,
     *                  adummy(1),temp,fmaxl
      logical ltem1,ltem2,reform,fbind,cdone,fdone,eqdone
      data ntwo/2/
c
c    The above data statement is to fool some compilers
c    that do not like a hard number as the index of w in
c    "call resign ..."
c
      nlin=nnineq-nineqn
      ii=1
      itry=1
      steps=1.d0
      ostep=steps
      fbind=.false.
      cdone=.false.
      fdone=.false.
      eqdone=.false.
      if(local) dlfeas=.false.
      ikeep=nlin-iskp
      prod1=0.1d0*grdftd
      tolfe=0.d0
      if(lqpsl) tolfe=tolfea
      if(iprint.ge.3.and.ipyes.eq.0)
     *  call sbout1(io,0,'directional deriv.  ',grdftd,adummy,1,2)

      w(1)=fM
 100  continue
        reform=.true.
        if(iprint.ge.3.and.ipyes.eq.0)
     *    write(io,9901) itry
        prod=prod1*steps
        if(.not.feasb.or.nobL.gt.1) prod=prod+tolfe
        do 200 i=1,nparam
          if(local)      xnew(i)=x(i)+steps*di(i)
          if(.not.local) xnew(i)=x(i)+steps*di(i)+d(i)*steps**2
 200    continue
        xisnew=.true.
        if(iprint.lt.3.or.ipyes.gt.0) goto 205
          call sbout1(io,0,'trial step          ',steps,adummy,1,2)
          call sbout1(io,nparam,'trial point         ',
     *                dummy,xnew,2,2)
 205    if(iskp.eq.0) goto 209
          ostep=steps
          do 207 i=ii,iskp
            ij=iskip(i)
            call constr(nparam,ij,xnew,g(ij))
            if(iprint.lt.3.or.ipyes.gt.0) goto 206
              if(i.eq.1) write(io,9900) ij,g(ij)
              if(i.ne.1) write(io,9902) ij,g(ij)
 206        if(g(ij).le.tolfe) goto 207
            ii=i
            goto 1120
 207      continue
          iskp=0
 209    if(nn.eq.0) goto 310
        if(.not.local.and.fbind) goto 315
 210    continue
        do 300 i=1,nn
          ncg=i
          ii=iact(i)
          ij=nnineq+neqn
          if(ii.le.nineqn.and.istore(ii).eq.1) goto 215
          if(ii.gt.nnineq.and.ii.le.ij.and.eqdone) goto 215
            temp=1.d0
            if(ii.gt.nnineq.and.ii.le.ij) temp=signeq(ii-nnineq)
            call constr(nparam,ii,xnew,g(ii))
            g(ii)=g(ii)*temp
            ncallg=ncallg+1
 215      if(iprint.lt.3.or.ipyes.gt.0) goto 220
            if(i.eq.1.and.ikeep.eq.nlin)
     *        write(io,9900) ii,g(ii)
            if(i.ne.1.or.ikeep.ne.nlin) write(io,9902) ii,g(ii)
 220      if(local.or.g(ii).le.tolfe) goto 230
            call shift(nn,ii,iact)
            goto 1110
 230      if(local.and.g(ii).gt.tolfe) goto 1500
 300    continue
 310    cdone=.true.
        eqdone=.true.
        if(local) dlfeas=.true.
 315    if(fdone) goto 410
        i = 0
        if(nob.gt.0) fmaxl=-bigbnd
 400    continue
          if(i.gt.nob) goto 405
          if(nob.ne.0.and.i.eq.0) i = 1
          ncf=i
          ii=iact(nn+i)
          if(feasb) then
            if(eqdone.or.neqn.eq.0) goto 317
              do 316 j=1,neqn
 316            call constr(nparam,nnineq+j,xnew,g(nnineq+j))
              ncallg=ncallg+neqn
 317        if(neqn.eq.0) goto 318
              if(eqdone)      job=20
              if(.not.eqdone) job=10
              call resign(nparam,neqn,psf,w(ntwo),penp,
     *                    g(nnineq+1),w(ntwo),signeq,job,10)
 318        if(i.eq.0.or.istore(nineqn+ii).eq.1) goto 320
              call obj(nparam,ii,xnew,f(ii))
              ncallf=ncallf+1
 320        if(i.eq.0) fii = 0.d0
            if(i.ne.0) fii = f(ii)
            if(nob.gt.0.and.(i.le.1.and.iprint.ge.3.and.ipyes.eq.0))
     *        write(io,9903) ii,fii-psf
            if(nob.eq.0.and.(i.le.1.and.iprint.ge.3.and.ipyes.eq.0))
     *        write(io,9904) fii-psf
            if(i.gt.1.and.iprint.ge.3.and.ipyes.eq.0)
     *        write(io,9902) ii,fii-psf
          else
            if(istore(ii).eq.1) goto 325
              call constr(nparam,indxob(ii),xnew,f(ii))
              ncallg=ncallg+1
 325        if(f(ii).gt.tolfe) reform=.false.
            if(i.eq.1.and.iprint.ge.3.and.ipyes.eq.0)
     *        write(io,9903) indxob(ii),f(ii)
            if(i.ne.1.and.iprint.ge.3.and.ipyes.eq.0)
     *        write(io,9902) indxob(ii),f(ii)
            fii=f(ii)
          endif
          fmaxl=dmax1(fmaxl,fii)
          if(nobL.ne.nob) fmaxl=dmax1(fmaxl,-fii)
          if(.not.feasb.and.reform) goto 401
          if(local) goto 340
          if((fii-psf).le.(fMp+prod)) goto 330
            fbind=.true.
            call shift(nob,ii,iact(nn+1))
          goto 1110
 330      if(nobL.eq.nob.or.(-fii-psf).le.(fMp+prod)) goto 401
            fbind=.true.
            call shift(nob,ii,iact(nn+1))
          goto 1110
 340      ltem1=(fii-psf).gt.(fMp+prod)
          ltem2=nobL.ne.nob.and.(-fii-psf).gt.(fMp+prod)
          if(ltem1.or.ltem2) goto 1500
 401      i = i + 1
        goto 400
 405    fbind=.false.
        fdone=.true.
        eqdone=.true.
        if(.not.cdone) goto 210
 410    if(ostep.eq.steps) mnm=ikeep+neq-neqn
        if(ostep.ne.steps) mnm=ncnstr-nn
        do 500 i=1,mnm
          ii=indxcn(i+nn)
          if(ikeep.ne.nlin.and.ostep.eq.steps) then
            if(i.le.ikeep) ii=iskip(nlin+2-i)
            if(i.gt.ikeep) ii=indxcn(nn+i-ikeep+nlin)
          endif
          call constr(nparam,ii,xnew,g(ii))
 500    continue
        scvneq=0.d0
        do 600 i=1,ncnstr
          if(i.gt.nnineq.and.i.le.(nnineq+neqn)) scvneq=scvneq-g(i)
 600      backup(i)=g(i)
        do 700 i=1,nob
 700      backup(i+ncnstr)=f(i)
        if(feasb.or..not.reform) goto 810
          do 800 i=1,nparam
 800        x(i)=xnew(i)
          nstop=0
          goto 1500
 810    if(local) ncg=ncnstr
        if(local) update=.true.
        fM=fmaxl
        fMp=fmaxl-psf
        fmax=fmaxl
        do 1000 i=1,nn
 1000     iact(i)=indxcn(i)
        do 1100 i=1,nob
 1100     iact(nn+i)=i
        goto 1500
c
 1110   cdone=.false.
        fdone=.false.
        eqdone=.false.
        reform=.false.
        if(lstype.eq.2) fbind=.false.
 1120   itry=itry+1
        if(steps.lt.1.d0) goto 1140
        do 1130 i=1,nob+nineqn
 1130     istore(i)=0
 1140   steps=steps*.5d0
        if(steps.lt.epsmac) goto 1150
      goto 100
c
 1150 info=4
      nstop=0
 1500 if(steps.lt.6.d-1) goto 9000
        do 1600 i=1,nob+nineqn
 1600     istore(i)=0
 9000 return
 9900 format(1x,t17,17htrial constraints,t37,i7,t45,e22.14)
 9901 format(1x,t17,12htrial number,t45,i22)
 9902 format(1x,t37,i7,t45,e22.14)
 9903 format(1x,t17,16htrial objectives,t37,i7,t45,e22.14)
 9904 format(1x,t17,25htrial penalized objective,t45,e22.14)
      end

      subroutine hesian(nparam,nob,nobL,nineqn,neq,neqn,nn,ncnstr,
     *                  nctotl,nfs,nstart,feasb,bigbnd,xnew,x,f,fmax,
     *                  fM,fMp,psf,gradf,grdpsf,penp,g,gm,gradg,indxob,
     *                  indxcn,cllamd,delta,eta,gamma,hess,hd,steps,
     *                  nrst,signeq,span,obj,constr,gradob,gradcn,
     *                  phess,psb,psmu,w,lenw,iw,leniw)
c
c     FFSQP : updating the Hessian matrix using BFGS
c             formula with Powell's modification
c
c     implicit real*8(a-h,o-z)
      integer nparam,nob,nobL,nineqn,neq,neqn,nn,ncnstr,nctotl,nfs,
     *        nstart,indxob(1),indxcn(1),nrst,lenw,leniw,iw(leniw)
c    *        nstart,indxob(nob),indxcn(1),nrst,lenw,leniw,iw(leniw)
      double  precision bigbnd,steps,psf,fmax,fM,fMp,
     *        xnew(nparam),x(nparam),f(1),gradf(nparam,1),
     *        grdpsf(nparam),penp(1),g(1),gm(1),
     *        gradg(nparam,1),cllamd(1),delta(nparam),
     *        eta(nparam),gamma(nparam),hess(nparam,nparam),hd(nparam),
     *        signeq(1),span(1),phess(1),psb(1),psmu(1),w(lenw)
c     double  precision bigbnd,steps,psf,fmax,fM,fMp,
c    *        xnew(nparam),x(nparam),f(nob),gradf(nparam,nob),
c    *        grdpsf(nparam),penp(neqn),g(ncnstr),gm(4*neqn),
c    *        gradg(nparam,ncnstr),cllamd(nctotl),delta(nparam),
c    *        eta(nparam),gamma(nparam),hess(nparam,nparam),hd(nparam),
c    *        signeq(neqn),span(1),phess(neq,neq),psb(neq),
c    *        psmu(neq),w(lenw)
      external obj,constr,gradob,gradcn
      logical feasb
c
      integer nnineq,M,ncallg,ncallf,mode,io,iprint,ipspan,ipyes,info,
     *        ipd,iter,nstop,initvl,lstype
      double  precision epsmac,rteps,udelta,valnom,objeps,objrep,gLgeps
      logical rolis1,d0is0,xisnew
      common  /fsqpp1/nnineq,M,ncallg,ncallf,mode,lstype,nstop,
     *        /fsqpp2/io,iprint,ipspan,ipyes,info,ipd,iter,initvl,
     *        /fsqpp3/epsmac,rteps,udelta,valnom,
     *        /fsqpp4/rolis1,d0is0,
     *        /fsqpus/objeps,objrep,gLgeps,xisnew
c
      integer ng,i,j,ifail,indexs,np,mnm,iout
      double  precision dhd,gammd,etad,scaprd,dummy,theta,signgj,psfnew
      logical done
c
      if(feasb.and.nstop.ne.0.and.neqn.eq.0) then
c
c       check of gLgeps is just after computing d0!
c
        if(dabs(w(1)-fmax).le.objeps) then
          nstop=0
        else if(dabs(w(1)-fmax).le.objrep*dabs(w(1))) then
          nstop=0
        endif
      endif
      if(nstop.eq.0) goto 810
c
      ipd=0
      done=.false.
      psfnew=0.d0
      call nullvc(nparam,delta)
      call nullvc(nparam,eta)
      if(nobL.gt.1) ng=2
      if(nobL.le.1) ng=1
c
 100  continue
        call nullvc(nparam,gamma)
        if(nobL.gt.1) call matrvc(nparam,nob,nparam,nob,gradf,
     *                     cllamd(nparam+ng+ncnstr),hd)
        if(.not.feasb) goto 120
        if(nineqn.eq.0) goto 110
        call matrvc(nparam,nineqn,nparam,nineqn,gradg,cllamd(nparam+ng),
     *              gamma)
 110    if(neqn.eq.0) goto 120
        call matrvc(nparam,neqn,nparam,neqn,gradg(1,nnineq+1),
     *              cllamd(nparam+nnineq+ng),eta)
 120    do 200 i=1,nparam
          if(nobL.gt.1) then
            if(done) psb(i)=hd(i)+cllamd(i)+gamma(i)
            gamma(i)=gamma(i)+hd(i)-grdpsf(i)+eta(i)
          else if(nobL.eq.1) then
            if(done) psb(i)=gradf(i,1)+cllamd(i)+gamma(i)
            gamma(i)=gamma(i)+gradf(i,1)-grdpsf(i)+eta(i)
          else if(nobL.eq.0) then
            if(done) psb(i)=cllamd(i)+gamma(i)
            gamma(i)=gamma(i)-grdpsf(i)+eta(i)
          endif
          if(.not.done) delta(i)=gamma(i)
 200    continue
        if(done) goto 410
        if(d0is0) goto 405
        if(nn.eq.0) goto 310
        do 300 i=1,nn
          if(feasb.and.i.gt.nineqn)     signgj=signeq(i-nineqn)
          if(.not.feasb.or.i.le.nineqn) signgj=1.d0
          valnom=g(indxcn(i))*signgj
          call gradcn(nparam,indxcn(i),xnew,gradg(1,indxcn(i)),constr)
 300    continue
        call resign(nparam,neqn,psf,grdpsf,penp,g(nnineq+1),
     *              gradg(1,nnineq+1),signeq,11,11)
 310    do 400 i=1,nob
          valnom=f(i)
          if(feasb) call gradob(nparam,i,xnew,gradf(1,i),obj)
          if(.not.feasb)
     *      call gradcn(nparam,indxob(i),xnew,gradf(1,i),constr)
 400    continue
 405    done=.true.
      goto 100
c
 410  if(d0is0) goto 910
      if(nrst.lt.(5*nparam).or.steps.gt.0.1d0) goto 420
        nrst=0
        call diagnl(nparam,1.d0,hess)
        goto 810
 420  nrst=nrst+1
      do 500 i=1,nparam
        gamma(i)=gamma(i)-delta(i)
        delta(i)=xnew(i)-x(i)
 500  continue
      call matrvc(nparam,nparam,nparam,nparam,hess,delta,hd)
      dhd=scaprd(nparam,delta,hd)
      if(sqrt(scaprd(nparam,delta,delta)).le.epsmac) then
        nstop=0
        info=8
        goto 9000
      endif
      gammd=scaprd(nparam,delta,gamma)
      if(gammd.ge.0.2d0*dhd) theta=1.d0
      if(gammd.lt.0.2d0*dhd) theta=.8d0*dhd/(dhd-gammd)
      do 600 i=1,nparam
 600    eta(i)=hd(i)*(1.d0-theta)+theta*gamma(i)
      etad=theta*gammd+(1.d0-theta)*dhd
      do 800  i=1,nparam
        do 700 j=i,nparam
          hess(i,j)=hess(i,j)-hd(i)*hd(j)/dhd+eta(i)*eta(j)/etad
 700    hess(j,i)=hess(i,j)
 800  continue
 810  do 900 i=1,nparam
 900    x(i)=xnew(i)
      xisnew=.true.
 910  if(nstop.eq.0) goto 9000
      if(neqn.eq.0.or..not.feasb) goto 1400
        iout=io
        i=neq-neqn
        if(i.eq.0) goto 940
        call matrvc(nparam,i,nparam,i,gradg(1,nnineq+neqn+1),
     *              cllamd(nparam+ng+nnineq+neqn),gamma)
        do 930 i=1,nparam
 930      psb(i)=psb(i)+gamma(i)
 940    i=nnineq-nineqn
        if(i.eq.0) goto 990
        call matrvc(nparam,i,nparam,i,gradg(1,nineqn+1),
     *              cllamd(nparam+ng+nineqn),gamma)
        do 950 i=1,nparam
 950      psb(i)=psb(i)+gamma(i)
 990    call estlam(nparam,neqn,ifail,iout,bigbnd,phess,delta,eta,gamma,
     *              gradg(1,nnineq+1),psb,hd,xnew,psmu,w,lenw,iw,leniw)
        do 1000 i=1,neqn
          if(ifail.ne.0.or.d0is0) then
            penp(i)=2.d0*penp(i)
          else if(ifail.eq.0) then
            etad=psmu(i)+penp(i)
            if(etad.ge.1.d0) goto 1000
            penp(i)=dmax1(1.0d0-psmu(i),2.0d0*penp(i))
          endif
          if(penp(i).gt.bigbnd) then
            nstop=0
            info=9
            goto 9000
          endif
 1000   continue
        call resign(nparam,neqn,psf,grdpsf,penp,g(nnineq+1),
     *              gradg(1,nnineq+1),signeq,20,12)
        fMp=fM-psf
 1400   if(nfs.eq.0) goto 1430
        nstart=nstart+1
        np=indexs(nstart,nfs)
        span(np)=fmax
        do 1410 i=1,neqn
 1410     gm((np-1)*neqn+i)=g(nnineq+i)
        if(neqn.ne.0) call resign(nparam,neqn,psfnew,grdpsf,penp,
     *                            gm(1),gradg,signeq,20,10)
        fM=span(1)
        fMp=span(1)-psfnew
        mnm=min0(nstart,nfs)
        do 1420 i=2,mnm
          if(neqn.ne.0) call resign(nparam,neqn,psfnew,grdpsf,penp,
     *                           gm((i-1)*neqn+1),gradg,signeq,20,10)
          fM=dmax1(fM,span(i))
          fMp=dmax1(fMp,span(i)-psfnew)
 1420   continue
 1430 if(iprint.lt.3.or.ipyes.gt.0) goto 9000
        do 1700 i=1,nob
          if(.not.feasb) goto 1600
            if(nob.gt.1) call sbout2(io,nparam,i,'gradf(j,',')',
     *                               gradf(1,i))
            if(nob.eq.1) call sbout1(io,nparam,'gradf(j)            ',
     *                               dummy,gradf(1,i),2,2)
          goto 1700
 1600       call sbout2(io,nparam,indxob(i),'gradg(j,',')',
     *                  gradf(1,i))
 1700   continue
        if(ncnstr.eq.0) goto 1900
        do 1800 i=1,ncnstr
 1800     call sbout2(io,nparam,indxcn(i),'gradg(j,',')',
     *                gradg(1,indxcn(i)))
        if(neqn.eq.0) goto 1900
        call sbout1(io,nparam,'grdpsf(j)           ',dummy,grdpsf,2,2)
        call sbout1(io,neqn,'P                   ',dummy,penp,2,2)
c       call sbout1(io,neqn,'psmu                ',dummy,psmu,2,2)
 1900   call sbout1(io,nparam,'multipliers  for  x ',dummy,cllamd,2,2)
        if(ncnstr.ne.0) call sbout1(io,ncnstr,'             for  g ',
     *                              dummy,cllamd(nparam+ng),2,2)
        if(nobL.gt.1) call sbout1(io,nob,'             for  f ',
     *                            dummy,cllamd(nparam+ng+ncnstr),2,2)
        do 2000 i=1,nparam
 2000     call sbout2(io,nparam,i,'hess (j,',')',hess(1,i))
 9000 return
      end
      subroutine grobfd(nparam,j,x,gradf,obj)
c
c     FFSQP : computation of gradients of objective
c             functions by forward finite differences
c
c     implicit real*8(a-h,o-z)
      integer nparam,j
      double  precision x(nparam),gradf(nparam)
      external obj
c
      integer io,iprint,ipspan,ipyes,info,ipd,idum,idum2
      double  precision epsmac,rteps,udelta,fj,objeps,objrep,gLgeps
      logical xisnew
      common  /fsqpp2/io,iprint,ipspan,ipyes,info,ipd,idum,idum2,
     *        /fsqpp3/epsmac,rteps,udelta,fj
      common  /fsqpus/objeps,objrep,gLgeps,xisnew
c
c     estimates the gradient of the objective function
c     by forward finite differences
c
      integer i
      double  precision xi,delta,dmax1
c
      do 10 i=1,nparam
        xi=x(i)
        delta=dmax1(udelta,rteps*dmax1(1.d0,dabs(xi)))
        if (xi.lt.0.d0) delta=-delta
        if (ipd.eq.1.or.j.ne.1.or.iprint.lt.3.or.ipyes.gt.0) goto 9
          if(i.eq.1) write(io,1001) delta
          if(i.ne.1) write(io,1002) delta
  9     x(i)=xi+delta
        xisnew=.true.
        call obj(nparam,j,x,gradf(i))
        gradf(i)=(gradf(i)-fj)/delta
        x(i)=xi
        xisnew=.true.
 10   continue
      return
 1001 format(1x,t17,8hdelta(i),t45,e22.14)
 1002 format(1x,t45,e22.14)
      end
c
      subroutine grcnfd(nparam,j,x,gradg,constr)
c
c     FFSQP : computation of gradients of constraint
c             functions by forward finite differences
c
c     implicit real*8(a-h,o-z)
      integer nparam,j
      double  precision x(nparam),gradg(nparam)
      external constr
c
      integer io,iprint,ipspan,ipyes,info,ipd,idum,idum2
      double  precision epsmac,rteps,udelta,gj,objeps,objrep,gLgeps
      logical xisnew
      common  /fsqpp2/io,iprint,ipspan,ipyes,info,ipd,idum,idum2,
     *        /fsqpp3/epsmac,rteps,udelta,gj
      common  /fsqpus/objeps,objrep,gLgeps,xisnew
c
c     estimate the gradient of the ith constraint
c     by forward finite differences
c
      integer i
      double  precision xi,delta,dmax1
c
      do 10 i=1,nparam
        xi=x(i)
        delta=dmax1(udelta,rteps*dmax1(1.d0,dabs(xi)))
        if (xi.lt.0.d0) delta=-delta
        if (j.ne.1.or.iprint.lt.3) goto 9
        if (ipspan.ge.10.and.ipyes.gt.0) goto 9
          if(i.eq.1) write(io,1001) delta
          if(i.ne.1) write(io,1002) delta
        ipd=1
  9     x(i)=xi+delta
        xisnew=.true.
        call constr(nparam,j,x,gradg(i))
        gradg(i)=(gradg(i)-gj)/delta
        x(i)=xi
        xisnew=.true.
 10   continue
      return
 1001 format(1x,t17,8hdelta(i),t45,e22.14)
 1002 format(1x,t45,e22.14)
      end
      subroutine out(miter,nparam,nob,nineqn,nn,neqn,ncnstr,x,g,f,
     *               fmax,fM,psf,steps,sktnom,d0norm,feasb)
c
c     FFSQP : output for different value of iprint
c
c     implicit real*8(a-h,o-z)
      integer miter,nparam,nob,nineqn,nn,neqn,ncnstr
      double  precision fmax,fM,steps,sktnom,d0norm,psf
      double  precision x(nparam),g(1),f(1)
c     double  precision x(nparam),g(ncnstr),f(nob)
      logical feasb
c
      integer nnineq,M,ncallg,ncallf,mode,io,iprint,ipspan,ipyes,
     *        info,idum1,iter,nstop,initvl,lstype
      common  /fsqpp1/nnineq,M,ncallg,ncallf,mode,lstype,nstop,
     *        /fsqpp2/io,iprint,ipspan,ipyes,info,idum1,iter,initvl
c
      integer i
      double precision SNECV,dummy,adummy(1)
c
      if(nstop.eq.0) ipyes=0
      if (iter.ge.miter.and.nstop.ne.0) then
        info=3
        nstop=0
      endif
 10   if(iprint.eq.0.or.ipyes.gt.0) then
        iter=iter+1
        goto 9000
      endif
      if(info.gt.0.and.info.lt.3) goto 120
      if(iprint.ne.1.or.nstop.eq.0) goto 20
        iter=iter+1
        if(initvl.eq.0) goto 9000
        if(feasb.and.nob.gt.0)
     *    call sbout1(io,nob,'objectives          ',dummy,f,2,1)
        if (mode.eq.1.and.iter.gt.1.and.feasb)
     *  call sbout1(io,0,'objective max4      ',fM,adummy,1,1)
        if(nob.gt.1) call sbout1(io,0,'objmax              ',
     *                           fmax,adummy,1,1)
        if(ncnstr.eq.0) write(io,9909)
        call sbout1(io,ncnstr,'constraints         ',dummy,g,2,1)
        if(ncnstr.ne.0) write(io,9909)
        goto 9000
 20   if(iprint.eq.1.and.nstop.eq.0) write(io,9900) iter
      if(iprint.ge.2.and.nstop.eq.0.and.ipspan.ge.10)
     *  write(io,9900) iter
      iter=iter+1
      if(initvl.eq.0)
     *  call sbout1(io,nparam,'x                   ',dummy,x,2,1)
      if(nob.gt.0)
     *  call sbout1(io,nob,'objectives          ',dummy,f,2,1)
      if (mode.eq.1.and.iter.gt.1)
     *  call sbout1(io,0,'objective max4      ',fM,adummy,1,1)
      if(nob.gt.1) call sbout1(io,0,'objmax              ',
     *                          fmax,adummy,1,1)
      if(ncnstr.eq.0) go to 110
      call sbout1(io,ncnstr,'constraints         ',dummy,g,2,1)
      if(.not.feasb) goto 110
        SNECV=0.d0
        do 100 i=nnineq+1,nnineq+nn-nineqn
          SNECV=SNECV+dabs(g(i))
 100    continue
        if(initvl.eq.0.and.(nn-nineqn).ne.0)
     *    call sbout1(io,0,'SNECV               ',SNECV,adummy,1,1)
 110  continue
      if(iter.le.1) write(io,9909)
      if(iter.le.1.and.ipspan.lt.10) write(io,9900) iter
      if(iter.le.1) goto 9000
      if(iprint.ge.2.and.initvl.eq.0)
     *  call sbout1(io,0,'step                ',steps,adummy,1,1)
      if(initvl.eq.0.and.(nstop.eq.0.or.info.ne.0.or.iprint.eq.2)) then
        call sbout1(io,0,'d0norm              ',d0norm,adummy,1,1)
        call sbout1(io,0,'ktnorm              ',sktnom,adummy,1,1)
      endif
      if(initvl.eq.0.and.feasb) write(io,9902) ncallf
      if(initvl.eq.0.and.(nn.ne.0.or..not.feasb)) write(io,9903) ncallg
      if(nstop.ne.0) write(io,9909)
      if(nstop.ne.0.and.iter.le.miter.and.ipspan.lt.10)
     *  write(io,9900) iter
 120  if(nstop.ne.0.or.iprint.eq.0) goto 9000
      write(io,9909)
      write(io,9901) info
      if(info.eq.0) write(io,9904)
      if(info.eq.0.and.sktnom.gt.0.1d0) write(io,9910)
      if(info.eq.3) write(io,9905)
      if(info.eq.4) write(io,9906)
      if(info.eq.5) write(io,9907)
      if(info.eq.6) write(io,9908)
      if(info.eq.8) write(io,9911)
      if(info.eq.9) write(io,9912)
      write(io,9909)
 9000 initvl=0
      if(ipspan.ge.10) ipyes=mod(iter,ipspan)
c      if(iter.le.miter) return
c        nstop=0
c        info=3
c        write(io,9905)
      return
 9900 format(1x,9hiteration,t22,i22)
 9901 format(1x,6hinform,t22,i22)
 9902 format(1x,6hncallf,t22,i22)
 9903 format(1x,6hncallg,t22,i22)
 9904 format(1x,'Normal termination: You have obtained a solution !!')
 9905 format(1x,'Warning : Maximum iterations have been reached ',
     *          'before obtaining a solution !!'/)
 9906 format(1x,'Error : Step size has been smaller than ',
     *          'the computed machine precision !!'/)
 9907 format(1x,'Error : Failure of the QP solver ',
     *          'in constructing d0 !!',
     *      /1x,'        A more robust QP solver may succeed.'/)
 9908 format(1x,'Error : Failure of the QP solver ',
     *          'in constructing d1 !!',
     *      /1x,'        A more robust QP solver may succeed.'/)
 9909 format(1x,/)
 9910 format(1x,'Warning: Norm of Kuhn-Tucker vector is large !!'/)
 9911 format(1x,'Error : The new iterate is numerically equivalent to ',
     *      /1x,'the current iterate, though the stopping criterion',
     *      /1x,'is not satisfied.'/)
 9912 format(1x,'Error : Could not satisfy nonlinear equality',
     *      /1x,'constraints - penalty parameter too large.'/)
      end
c=== subroutines used in FFSQP ====================================c
c                                                                  c
c  diagnl  error   estlam  fool    indexs  lfuscp  matrcp  matrvc  c
c  nullvc  resign  sbout1  sbout2  scaprd  shift   slope   small   c
c                                                                  c
c==================================================================c
c
      subroutine diagnl(nrowa,diag,a)
c     implicit real*8(a-h,o-z)
      integer nrowa,i,j
      double  precision a(nrowa,1),diag
c     double  precision a(nrowa,nrowa),diag
c
c     set a=diag*I, the diagonal matrix
c
      do 200 i=1,nrowa
        do 100 j=i,nrowa
          a(i,j)=0.d0
 100      a(j,i)=0.d0
 200    a(i,i)=diag
      return
      end
c
      subroutine error(string,inform,io)
c     implicit real*8 (a-h,o-z)
      integer inform,io
      character*40 string
c
      write(io,9900) string
 9900 format(1x,a40)
      inform=7
      return
      end
c
      subroutine estlam(nparam,neqn,ifail,iout,bigbnd,hess,cvec,a,b,
     *                  gradh,psb,bl,bu,x,w,lenw,iw,leniw)
      integer nparam,neqn,ifail,iout,lenw,leniw,iw(leniw)
      double precision bigbnd,hess(neqn,1),cvec(1),a(1),b(1),
     *                 gradh(nparam,1),psb(1),bl(1),bu(1),
     *                 x(1),w(lenw)
c     double precision bigbnd,hess(neqn,neqn),cvec(neqn),a(1),b(1),
c    *                 gradh(nparam,neqn),psb(nparam),bl(1),bu(1),
c    *                 x(neqn),w(lenw)
c
c     compute an estimate of multipliers for updating penalty parameter
c
      integer i,j
      double precision scaprd
c
      do 200 i=1,neqn
        bl(i)=-bigbnd
        bu(i)=bigbnd
        cvec(i)=scaprd(nparam,gradh(1,i),psb)
        x(i)=0.d0
        do 100 j=i,neqn
          hess(i,j)=scaprd(nparam,gradh(1,i),gradh(1,j))
 100      hess(j,i)=hess(i,j)
 200  continue
      iw(1)=1
      call ql0001(0,0,1,neqn,neqn,2*neqn,hess,cvec,a,b,bl,bu,x,w,
     c            iout,ifail,0,w(2),lenw-1,iw,leniw)
c
      return
      end
c
      subroutine fool(x,y,z)
      double precision x,y,z
c
      z=x*y+y
      return
      end
c
      double precision function lfuscp(val,thrshd)
c     implicit real*8(a-h,o-z)
      double precision val,thrshd
c
      if(dabs(val).le.thrshd) lfuscp=0
      if(dabs(val).gt.thrshd) lfuscp=1
      return
      end
c
      integer function indexs(i,nfs)
c     implicit real*8(a-h,o-z)
      integer i,nfs,mm
c
c     find the residue of i with respect to nfs
c
      mm=i
      if(mm.le.nfs) goto 120
 110  mm=mm-nfs
      if(mm.gt.nfs) goto 110
 120  indexs=mm
      return
      end
c
      subroutine matrcp(ndima,a,ndimb,b)
c     implicit real*8(a-h,o-z)
      integer ndima,ndimb,i,j
      double  precision a(ndima,1),b(ndimb,1)
c     double  precision a(ndima,ndima),b(ndimb,ndimb)
c
      do 100 i=1,ndima
        do 100 j=1,ndima
 100      b(i,j)=a(i,j)
      if(ndimb.le.ndima) goto 9000
        do 200 i=1,ndimb
          b(ndimb,i)=0.d0
 200      b(i,ndimb)=0.d0
 9000 return
      end
c
      subroutine matrvc(l,n,la,na,a,x,y)
c     implicit real*8(a-h,o-z)
      integer l,n,la,na,i,j
      double  precision a(l,n),x(n),y(l),yi
c     double  precision a(l,1),x(1),y(1),yi
c
c     computation of y=ax
c
      do 200 i=1,la
        yi=0.d0
        do 100 j=1,na
 100      yi=yi+a(i,j)*x(j)
 200      y(i)=yi
      return
      end
c
      subroutine nullvc(nparam,x)
c     implicit real*8(a-h,o-z)
      integer nparam,i
      double  precision x(nparam)
c
c     set x=0
c
      do 100 i=1,nparam
 100    x(i)=0.d0
      return
      end
c
      subroutine resign(n,neqn,psf,grdpsf,penp,g,gradg,signeq,job1,job2)
      integer i,j,job1,job2,n,neqn
      double precision psf,grdpsf(1),penp(1),g(1),gradg(n,1),
     *                 signeq(1)
c     double precision psf,grdpsf(n),penp(neqn),g(neqn),gradg(n,neqn),
c    *                 signeq(neqn)
c
c     job1=10: g*signeq, job1=11: gradg*signeq, job1=12: job1=10&11
c     job1=20: do not change sign
c     job2=10: psf,      job2=11: grdpsf,       job2=12: job2=10&11
c     job2=20: do not compute psf or grdpsf
c
      if(job2.eq.10.or.job2.eq.12) psf=0.d0
      do 100 i=1,neqn
        if(job1.eq.10.or.job1.eq.12) g(i)=signeq(i)*g(i)
        if(job2.eq.10.or.job2.eq.12) psf=psf+g(i)*penp(i)
        if(job1.eq.10.or.job1.eq.20) goto 100
          do 50 j=1,n
            gradg(j,i)=gradg(j,i)*signeq(i)
  50      continue
 100  continue
      if(job2.eq.10.or.job2.eq.20) goto 9000
      call nullvc(n,grdpsf)
      do 120 i=1,n
        do 110 j=1,neqn
 110      grdpsf(i)=grdpsf(i)+gradg(i,j)*penp(j)
 120  continue
c
 9000 return
      end
c
      subroutine sbout1(io,n,s1,z,z1,job,level)
c     implicit real*8(a-h,o-z)
      integer io,n,job,level,j
      double precision z,z1(1)
      character*20 s1
c
      if (job.eq.2) goto 10
      if (level.eq.1)write(io,9900) s1,z
      if (level.eq.2)write(io,9901) s1,z
      return
 10   if(n.eq.0) goto 101
      if (level.eq.1)write(io,9900) s1,z1(1)
      if (level.eq.2)write(io,9901) s1,z1(1)
      do 100 j=2,n
        if (level.eq.1) write(io,9902) z1(j)
        if (level.eq.2) write(io,9903) z1(j)
 100  continue
 101  return
 9900 format(1x,a20,e22.14)
 9901 format(1x,t17,a20,t45,e22.14)
 9902 format(1x,t22,e22.14)
 9903 format(1x,t45,e22.14)
      end
c
      subroutine sbout2(io,n,i,s1,s2,z)
c     implicit real*8(a-h,o-z)
      integer io,n,i,j
      double precision z(n)
      character*8 s1
      character*1 s2
c
      write(io,9900) s1,i,s2,z(1)
      do 100 j=2,n
 100    write(io,9901) z(j)
      return
 9900 format(1x,t17,a8,i5,a1,t45,e22.14)
 9901 format(1x,t45,e22.14)
      end
c
      double  precision function scaprd(n,x,y)
c     implicit real*8(a-h,o-z)
      integer n,i
      double  precision x(1),y(1),z
c     double  precision x(n),y(n),z
c
c     compute z=x'y
c
      z=0.d0
      do 100 i=1,n
        z=x(i)*y(i)+z
 100  continue
      scaprd=z
      return
      end
c
      subroutine shift(n,ii,iact)
      integer n,ii,iact(1),j,k
c
      if(ii.eq.iact(1)) return
      do 200 j=1,n
        if(ii.ne.iact(j)) goto 200
        do 100 k=j,2,-1
 100      iact(k)=iact(k-1)
        goto 210
 200  continue
 210  if(n.ne.0) iact(1)=ii
      return
      end
c
      double precision function slope(nob,nobL,neqn,nparam,feasb,
     *                         f,gradf,grdpsf,x,y,fmax,theta,job)
c     implicit real*8(a-h,o-z)
      integer nob,nobL,neqn,nparam,job,i
      double  precision fmax,theta,slope1,dmax1,dmin1,rhs,rhog,
     *        grdftx,grdfty,diff,scaprd,grpstx,grpsty
      double  precision f(nob),gradf(nparam,nob),grdpsf(nparam),
     *        x(nparam),y(nparam)
c     double  precision f(1),gradf(nparam,1),grdpsf(nparam),
c    *        x(nparam),y(nparam)
      logical feasb
c
      double  precision bigbnd,dummy
      common  /fsqpq1/bigbnd,dummy
c
c     job=0 : compute the generalized gradient of the minimax
c     job=1 : compute rhog in mode = 1
c
      slope=-bigbnd
      if(feasb.and.nob.eq.0) slope=0.d0
      if(neqn.eq.0.or..not.feasb) then
        grpstx=0.d0
        grpsty=0.d0
      else
        grpstx=scaprd(nparam,grdpsf,x)
        grpsty=scaprd(nparam,grdpsf,y)
      endif
      do 100 i=1,nob
        slope1=f(i)+scaprd(nparam,gradf(1,i),x)
        slope=dmax1(slope,slope1)
        if(nobL.ne.nob) slope=dmax1(slope,-slope1)
 100  continue
      slope=slope-fmax-grpstx
      if (job.eq.0) goto 9000
      rhs=theta*slope+fmax
      rhog=1.d0
      do 200 i=1,nob
        grdftx=scaprd(nparam,gradf(1,i),x)-grpstx
        grdfty=scaprd(nparam,gradf(1,i),y)-grpsty
        diff=grdfty-grdftx
        if (diff.le.0.d0) goto 200
        rhog=dmin1(rhog,(rhs-f(i)-grdftx)/diff)
        if(nobL.ne.nob) rhog=dmin1(rhog,-(rhs+f(i)+grdftx)/diff)
 200  continue
      slope=rhog
 9000 return
      end
c
      double precision function small()
c     implicit real*8(a-h,o-z)
      double precision one, two, z
c
      one=1.d0
      two=2.d0
      small=one
10    small=small/two
      call fool(small,one,z)
      if(z.gt.one) goto 10
      small=small*two*two
c
c The simpler sequence commented out below fails on some machines that use
c extra-length registers for internal computation.  This was pointed out
c to us by Roque Donizete de Oliveira (Michigan) who suggested to sequence
c used now.
c
c     small=1.d0
c100  if ((small+1.d0).eq.1.d0) goto 110
c     small=small/2.d0
c     goto 100
c110  small=small*4.d0
      return
      end
c
      block data
      double  precision objeps,objrep,gLgeps
      logical xisnew
      common  /fsqpus/objeps,objrep,gLgeps,xisnew
c
      data objeps,objrep,gLgeps/-1.d0,-1.d0,-1.d0/
      end
