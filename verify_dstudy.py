import numpy as np
from math import pi, exp
from scipy.stats import binom, logistic
plogis = logistic.cdf

def dstudy_bernoulli_draws(alpha, sd_obj, sd_facet, n_grid, K=500, K_facet=500, info=False, seed=None):
    rng = np.random.default_rng(seed)
    alpha=np.atleast_1d(alpha); sd_obj=np.atleast_1d(sd_obj); sd_facet=np.atleast_1d(sd_facet)
    n_grid=np.asarray(n_grid,float); S=len(alpha); nG=len(n_grid)
    lrv=pi**2/3
    link=np.full((S,nG),np.nan); resp=np.full((S,nG),np.nan)
    inf=np.full((S,nG),np.nan) if info else None
    for s in range(S):
        a,so,sf=alpha[s],sd_obj[s],sd_facet[s]
        link[s]=so**2/(so**2+(sf**2+lrv)/n_grid)
        u=rng.normal(0,so,K)
        if sf>0:
            v=rng.normal(0,sf,(K,K_facet)); p_cond=plogis(a+u[:,None]+v).mean(1)
        else:
            p_cond=plogis(a+u)
        var_pi=p_cond.var(ddof=1); E=np.mean(p_cond*(1-p_cond))
        resp[s]=var_pi/(var_pi+E/n_grid)
        if info:
            for g,n in enumerate(n_grid.astype(int)):
                supp=np.arange(n+1)
                pmf=binom.pmf(supp[None,:],n,p_cond[:,None])
                Hc=np.mean(-np.sum(pmf*np.log(pmf+1e-300),1))
                pm=pmf.mean(0); Hm=-np.sum(pm*np.log(pm+1e-300))
                I=max(Hm-Hc,0); inf[s,g]=1-exp(-2*I)
    out={'link':link,'response':resp}
    if info: out['info']=inf
    return out

def find_n(mat,n_grid,th):
    med=np.median(mat,0); idx=np.where(med>=th)[0]
    return None if len(idx)==0 else n_grid[idx[0]]

fails=[]
def check(cond,name):
    print(("PASS " if cond else "FAIL ")+name); 
    if not cond: fails.append(name)

# T1 shape
o=dstudy_bernoulli_draws([0,.5],[1,1.2],[.5,.3],[1,2,5,10,20],K=200,K_facet=100,seed=1)
check(o['link'].shape==(2,5) and o['response'].shape==(2,5),"shape")
check(((o['link']>=0)&(o['link']<=1)).all() and ((o['response']>=0)&(o['response']<=1)).all(),"bounded [0,1]")

# T2 link closed form exact
so,sf=1.3,.6; ng=np.arange(1,26)
o=dstudy_bernoulli_draws(0,so,sf,ng,K=50,K_facet=50,seed=2)
exp_=so**2/(so**2+(sf**2+pi**2/3)/ng)
check(np.allclose(o['link'][0],exp_,atol=1e-12),"link closed form exact")

# T3 monotone, ->1
ng=[1,2,4,8,16,32,64,128,256]
o=dstudy_bernoulli_draws(.3,1,.4,ng,K=2000,K_facet=200,seed=3)
check((np.diff(o['link'][0])>=0).all(),"link monotone")
check((np.diff(o['response'][0])>=0).all(),"response monotone")
check(o['link'][0,-1]>.95 and o['response'][0,-1]>.95,"both ->1 at n=256")

# T4 n=1 response == binomial engagement ICC (same draws)
rng=np.random.default_rng(4); so,a,K=1.1,-.2,20000
u=rng.normal(0,so,K); p=plogis(a+u); ref=p.var(ddof=1)/(p.var(ddof=1)+np.mean(p*(1-p)))
o=dstudy_bernoulli_draws(a,so,0,[1],K=K,K_facet=1,seed=4)
check(abs(o['response'][0,0]-ref)<1e-8,f"n=1 response == engagement ICC ({o['response'][0,0]:.6f} vs {ref:.6f})")

# T5 response <= link
ok=True
for so in [.5,1,2]:
    for sf in [0,.5]:
        o=dstudy_bernoulli_draws(0,so,sf,[1,5,20],K=3000,K_facet=200,seed=5)
        ok&=(o['response'][0]<=o['link'][0]+1e-6).all()
check(ok,"response <= link at all n, all settings")

# T6 info curve
o=dstudy_bernoulli_draws(0,1,.3,[1,2,5,10,25],K=800,K_facet=100,info=True,seed=6)
check(((o['info']>=0)&(o['info']<=1)).all(),"info bounded")
check((np.diff(o['info'][0])>=-.01).all(),"info monotone (MC tol)")
check((o['info'][0]<=o['link'][0]+.02).all(),"info <= link (Thm 5)")

# T7 required n
ng=np.arange(1,41)
o=dstudy_bernoulli_draws(0,1,.5,ng,K=500,K_facet=100,seed=7)
req={th:(find_n(o['link'],ng,th),find_n(o['response'],ng,th)) for th in [.7,.8,.9]}
L=[req[t][0] for t in [.7,.8,.9]]; R=[req[t][1] for t in [.7,.8,.9]]
check(all(np.diff([x for x in L if x])>=0),"required n link non-decreasing in threshold")
check(all(np.diff([x for x in R if x])>=0),"required n response non-decreasing")
check(all((r>=l) for l,r in zip(L,R) if l and r),"response needs >= link")
print("required n:",req)

print("\n"+("ALL PASS" if not fails else f"FAILURES: {fails}"))
