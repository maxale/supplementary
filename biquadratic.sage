'''
Supplementary SageMath code for the paper:
* M. A. Alekseyev and Sz. Tengely. "On Integral Points on Biquadratic Curves and Near-Multiples of Squares in Lucas Sequences". Journal of Integer Sequences 17:6 (2014), Article 14.6.6. 
  URL: https://cs.uwaterloo.ca/journals/JIS/VOL17/Alekseyev/alekseyev3.html
  arXiv:1306.0883 [math.NT] https://doi.org/10.48550/arXiv.1306.0883

Originally available at http://shrek.unideb.hu/~tengely/biquadratic.sage
'''

def alldivisors(n):
    d1=n.divisors()
    return d1+[-d for d in d1]

def QuarticEQ(A):
    if A[0]<0:
        Rx.<U>=QQ[]
        rroots=(A[0]*U^4+A[1]*U^2+A[2]).real_roots()
        if len(rroots)==0:
            return []
        else:
            Umin=round(min(rroots))-1
            Umax=round(max(rroots))+1
            xposs=[k for k in [Umin..Umax] if (A[0]*k^4+A[1]*k^2+A[2]).is_square()]
            return [[k,(A[0]*k^4+A[1]*k^2+A[2]).isqrt()] for k in xposs]+[[k,-(A[0]*k^4+A[1]*k^2+A[2]).isqrt()] for k in xposs]
    if A[0].is_square():
        a0=A[0].isqrt()
        Rx.<U>=QQ[]
        if A[1]%(2*a0)==0:
            rroots=(2*a0*U^2+A[2]-1+A[1]//a0).real_roots()+(-2*a0*U^2+A[2]-1-A[1]//a0).real_roots()
            if len(rroots)==0:
                return []
            else:
                Umin=round(min(rroots))-1
                Umax=round(max(rroots))+1
                xposs=[k for k in [Umin..Umax] if (A[0]*k^4+A[1]*k^2+A[2]).is_square()]
                return [[k,(A[0]*k^4+A[1]*k^2+A[2]).isqrt()] for k in xposs]+[[k,-(A[0]*k^4+A[1]*k^2+A[2]).isqrt()] for k in xposs]
        else:
            poly1=expand(A[0]*U^4+A[1]*U^2+A[2]-(a0*U^2+floor(A[1]/(2*a0)))^2)
            poly2=expand(A[0]*U^4+A[1]*U^2+A[2]-(a0*U^2+ceil(A[1]/(2*a0)))^2)
            rroots=poly1.real_roots()+poly2.real_roots()
            if len(rroots)==0:
                return []
            else:
                Umin=round(min(rroots))-1
                Umax=round(max(rroots))+1
                xposs=[k for k in [Umin..Umax] if (A[0]*k^4+A[1]*k^2+A[2]).is_square()]
                return [[k,(A[0]*k^4+A[1]*k^2+A[2]).isqrt()] for k in xposs]+[[k,-(A[0]*k^4+A[1]*k^2+A[2]).isqrt()] for k in xposs]
    SFA1=squarefree_part(4*A[0]*A[2]-A[1]^2)
    SFA2=((4*A[0]*A[2]-A[1]^2)//SFA1).isqrt()
    SFB1=squarefree_part(-A[0])
    SFB2=((-A[0])//SFB1).isqrt()
    C=Conic(QQ,[1,SFA1,SFB1])
    print(C)
    if C.has_rational_point():
        if A[2].is_square():
            Sol=[[0,A[2].isqrt()],[0,-A[2].isqrt()]]
        else:
            Sol=[]
        Cpar=C.parametrization(morphism = False)
        Rx.<U>=QQ[]
        Pxy.<x,y>=QQ[]
        Mpar=matrix(ZZ,3,3,[Cpar[0][0].monomial_coefficient(x^2),Cpar[0][0].monomial_coefficient(x*y),Cpar[0][0].monomial_coefficient(y^2),
        Cpar[0][1].monomial_coefficient(x^2),Cpar[0][1].monomial_coefficient(x*y),Cpar[0][1].monomial_coefficient(y^2),
        Cpar[0][2].monomial_coefficient(x^2),Cpar[0][2].monomial_coefficient(x*y),Cpar[0][2].monomial_coefficient(y^2)])
        ldiv=alldivisors(max(Mpar.elementary_divisors()))
        Pxyz.<X, Y, Z> = QQ[]
        cf1=Cpar[0][0].monomial_coefficient(x^2)
        cf2=Cpar[0][0].monomial_coefficient(x*y)
        cf3=Cpar[0][0].monomial_coefficient(y^2)
        alfa1=Cpar[0][1].monomial_coefficient(x^2)
        alfa2=Cpar[0][1].monomial_coefficient(x*y)
        alfa3=Cpar[0][1].monomial_coefficient(y^2)
        F1=cf1*X^2+cf2*X*Y+cf3*Y^2
        F2=alfa1*X^2+alfa2*X*Y+alfa3*Y^2
        alfa=squarefree_part(4*alfa1*alfa3-alfa2^2)
        Kfield.<beta>=NumberField(U^2+alfa)
        TSlist=[]
        for t in ldiv:
            for s in divisors(SFA2*t):
                Ctemp=Conic(QQ,SFA2*F1-A[1]*F2-((2*A[0]*SFA2*t)//s)*Z^2)
                if Ctemp.has_rational_point() and (4*alfa1*SFA2*t//s).is_norm(Kfield):
                    TSlist.append([s,t])
        TSset=Set([(SFA2*w[1])//w[0] for w in TSlist])
        rpxy=SFA2*Cpar[0][0]-A[1]*Cpar[0][1]
        a1=rpxy.monomial_coefficient(x^2)
        a2=rpxy.monomial_coefficient(x*y)
        a3=rpxy.monomial_coefficient(y^2)
        print(TSlist)
        print(TSset)
        print(rpxy)
        for st in TSset:
            if not a1==0:
                cV=squarefree_part(4*a1*a3-a2^2)
                cX=squarefree_part(8*A[0]*st*a1)
                mV=Integer((4*a1*a3-a2^2)//cV).isqrt()
                mX=Integer((8*A[0]*st*a1)//cX).isqrt()
                C1=Conic(QQ,X^2+cV*Y^2-cX*Z^2)
                print(C1)
                if C1.has_rational_point():
                    C1par=C1.parametrization(morphism = False)
                    M1par=matrix(ZZ,3,3,[C1par[0][0].monomial_coefficient(x^2),C1par[0][0].monomial_coefficient(x*y),C1par[0][0].monomial_coefficient(y^2),
                    C1par[0][1].monomial_coefficient(x^2),C1par[0][1].monomial_coefficient(x*y),C1par[0][1].monomial_coefficient(y^2),
                    C1par[0][2].monomial_coefficient(x^2),C1par[0][2].monomial_coefficient(x*y),C1par[0][2].monomial_coefficient(y^2)])
                    l1div=alldivisors(max(M1par.elementary_divisors()))
                    T1=C1par[0][0]
                    T2=C1par[0][1]
                    T3=C1par[0][2]
                    ThuePoly=alfa1*mV^2*T1^2+(2*alfa2*a1*mV-2*alfa1*a2*mV)*T1*T2+(alfa1*a2^2-2*a1*a2*alfa2+4*a1^2*alfa3)*T2^2
                    T_init=gp.thueinit(ThuePoly.subs({y:1}),1)
                    for t1 in l1div:
                        for t2 in Integer(2*a1*mV).divisors():
                            T_solve=gp.thue(T_init,(4*st*a1^2*t1^2*mV^2)//(t2^2))
                            if len(T_solve)>0:
                                for w in T_solve:
                                    x_possible=(T3.subs(x=Integer(w[1]),y=Integer(w[2]))*t2/(mX*t1)).constant_coefficient()
                                    if x_possible.is_integral():
                                        xint=Integer(x_possible)
                                        if (A[0]*xint^4+A[1]*xint^2+A[2]).is_square():
                                            xy=[xint,(A[0]*xint^4+A[1]*xint^2+A[2]).isqrt()]
                                            xym=[xint,-(A[0]*xint^4+A[1]*xint^2+A[2]).isqrt()]
                                        if not xy in Sol:
                                            Sol.append(xy)
                                        if not xym in Sol:
                                            Sol.append(xym)
            elif not a3==0:
                cU=squarefree_part(4*a1*a3-a2^2)
                cX=squarefree_part(8*A[0]*st*a3)
                mU=Integer((4*a1*a3-a2^2)//cU).isqrt()
                mX=Integer((8*A[0]*st*a3)//cX).isqrt()
                C1=Conic(QQ,X^2+cU*Y^2-cX*Z^2)
                if C1.has_rational_point():
                    C1par=C1.parametrization(morphism = False)
                    M1par=matrix(ZZ,3,3,[C1par[0][0].monomial_coefficient(x^2),C1par[0][0].monomial_coefficient(x*y),C1par[0][0].monomial_coefficient(y^2),
                    C1par[0][1].monomial_coefficient(x^2),C1par[0][1].monomial_coefficient(x*y),C1par[0][1].monomial_coefficient(y^2),
                    C1par[0][2].monomial_coefficient(x^2),C1par[0][2].monomial_coefficient(x*y),C1par[0][2].monomial_coefficient(y^2)])
                    l1div=alldivisors(max(M1par.elementary_divisors()))
                    T1=C1par[0][0]
                    T2=C1par[0][1]
                    T3=C1par[0][2]
                    ThuePoly=alfa3*mU^2*T1^2+(2*alfa2*a3*mU-2*alfa3*a2*mU)*T1*T2+(alfa3*a2^2-2*a3*a2*alfa2+4*a3^2*alfa1)*T2^2
                    T_init=gp.thueinit(ThuePoly.subs({y:1}),1)
                    for t1 in l1div:
                        for t2 in Integer(2*a3*mU).divisors():
                            T_solve=gp.thue(T_init,(4*st*a3^2*t1^2*mU^2)//(t2^2))
                            if len(T_solve)>0:
                                for w in T_solve:
                                    x_possible=(T3.subs(x=Integer(w[1]),y=Integer(w[2]))*t2/(mX*t1)).constant_coefficient()
                                    if x_possible.is_integral():
                                        xint=Integer(x_possible)
                                        if (A[0]*xint^4+A[1]*xint^2+A[2]).is_square():
                                            xy=[xint,(A[0]*xint^4+A[1]*xint^2+A[2]).isqrt()]
                                            xym=[xint,-(A[0]*xint^4+A[1]*xint^2+A[2]).isqrt()]
                                        if not xy in Sol:
                                            Sol.append(xy)
                                        if not xym in Sol:
                                            Sol.append(xym)
            else:
                sfpart=Integer(2*A[0]*st*a2).squarefree_part()
                for t1 in alldivisors(sfpart):
                    t2=sfpart//t1
                    ThuePoly=alfa1*t1^2*x^4+alfa2*t1*t2*x^2+alfa3*t2^2
                    T_init=gp.thueinit(ThuePoly)
                    T_solve=gp.thue(T_init,st)
                    if len(T_solve)>0:
                        for w in T_solve:
                            if (not w[1]==0) and (not w[2]==0):
                                x2_possible=a2*t1*t2*Integer(w[1])^2*Integer(w[2])^2/(2*A[0]*st)
                                if x2_possible.is_integral() and x2_possible.is_square():
                                    xint=Integer(x2_possible).isqrt()
                                    if (A[0]*xint^4+A[1]*xint^2+A[2]).is_square():
                                        xy1=[xint,(A[0]*xint^4+A[1]*xint^2+A[2]).isqrt()]
                                        xy1m=[xint,-(A[0]*xint^4+A[1]*xint^2+A[2]).isqrt()]
                                        xy2=[-xint,(A[0]*xint^4+A[1]*xint^2+A[2]).isqrt()]
                                        xy2m=[-xint,-(A[0]*xint^4+A[1]*xint^2+A[2]).isqrt()]
                                    if not xy1 in Sol:
                                        Sol.append(xy1)
                                    if not xy1m in Sol:
                                        Sol.append(xy1m)
                                    if not xy2 in Sol:
                                        Sol.append(xy2)
                                    if not xy2m in Sol:
                                        Sol.append(xy2m)
        return Sol
    else:
        X=ZZ['X'].0
        fpoly=A[0]*X^4+A[1]*X^2+A[2]
        rootsfpoly=fpoly.roots(multiplicities=False)
        if len(rootsfpoly)>0:
            return [[k,0] for k in rootsfpoly]
        else:
            return []

def biquadratic(A):
    if A[1]^2-4*A[0]*A[2]==0 or A[0]==0 or A[2]==0:
        return -1
    else:
        E=EllipticCurve([0,A[1],0,A[0]*A[2],0])
        Egens=E.gens(proof=False)
        if E.gens_certain():
            Sols=[]
            Eip=E.integral_points()
            for a in alldivisors(A[0]*A[2]):
                sols=[]
                c=(A[0]*A[2])//a
                sols.append([a,A[1],c])
                filterSol1=[Epoint for Epoint in Eip if Epoint[0] % a==0]
                filterSol2=[Epoint for Epoint in filterSol1 if (Epoint[0]//a).is_square()]
                for p1 in filterSol2:
                    if p1[0]==0:
                        if c.is_square():
                            sols.append([0,c.isqrt()])
                            sols.append([0,-1*(c.isqrt())])
                    else:
                        x1=Integer(p1[0]//a).isqrt()
                        y1=p1[1]//(a*x1)
                        sols.append([x1,y1])
                        sols.append([x1,-y1])
                        sols.append([-x1,y1])
                        sols.append([-x1,-y1])
                Sols.append(sols)
            return Sols
        else:
            Sols=[]
            for a in alldivisors(A[0]*A[2]):
                sols=[]
                c=(A[0]*A[2])//a
                sols.append([a,A[1],c])
                sols.append(["Thue"])
                xysols=QuarticEQ([a,A[1],c])
                if len(xysols)>0:
                    for xy in xysols:
                        sols.append(xy)
                else:
                    sols.append([])
                Sols.append(sols)
            return Sols
