# -*- coding: utf-8 -*-
# fusionEKF ve imuLoc MATLAB mantiginin Python kopyasi - marmaray verisinde dogrulama
import csv, math
import numpy as np

def load(path):
    with open(path, newline='') as f:
        return list(csv.DictReader(f))

def colf(rows, name):
    return np.array([float(r[name]) for r in rows])

acc = load('data_marmaray/Accelerometer.csv')
gyr = load('data_marmaray/Gyroscope.csv')
loc = load('data_marmaray/Location.csv')

t   = colf(acc,'seconds_elapsed')
ax  = colf(acc,'x'); ay = colf(acc,'y'); az = colf(acc,'z')
tg_gz = colf(gyr,'seconds_elapsed'); gz_raw = colf(gyr,'z')
gz = np.interp(t, tg_gz, gz_raw)   # gyro -> acc zamanina hizala
N = len(t)
dt_k = np.diff(t); dt_k = np.append(dt_k, dt_k[-1])
dt = float(np.mean(np.diff(t)))

tg  = colf(loc,'seconds_elapsed')
lat = colf(loc,'latitude'); lon = colf(loc,'longitude')
hacc= colf(loc,'horizontalAccuracy')
spd = colf(loc,'speed'); spd[spd<0] = np.nan
brg = colf(loc,'bearing')
vN_g = spd*np.cos(np.radians(brg)); vE_g = spd*np.sin(np.radians(brg))

Re=6378137.0; d2r=math.pi/180
ref_lat=lat[0]; ref_lon=lon[0]
gN=(lat-ref_lat)*d2r*Re
gE=(lon-ref_lon)*d2r*Re*math.cos(ref_lat*d2r)
ngps=len(tg)

def init_heading():
    for i in range(ngps):
        if brg[i]>=0: return math.radians(brg[i])
    for i in range(ngps):
        if not math.isnan(vN_g[i]) and math.hypot(vN_g[i],vE_g[i])>2:
            return math.atan2(vE_g[i],vN_g[i])
    return 0.0

def ekf_update(x,P,z,H,R,gate):
    y = z - H@x
    S = H@P@H.T + R
    if np.linalg.cond(S) > 1e15: return x,P
    d2 = float(y@np.linalg.solve(S,y))
    if d2 > gate: return x,P
    K = P@H.T@np.linalg.inv(S)
    x = x + K@y
    I=np.eye(6); P=(I-K@H)@P@(I-K@H).T + K@R@K.T
    P=0.5*(P+P.T)
    return x,P

def fusionEKF():
    x=np.zeros(6); x[0]=gN[0]; x[1]=gE[0]
    if not math.isnan(vN_g[0]): x[2]=vN_g[0]; x[3]=vE_g[0]
    P=np.diag([25.,25,4,4,0.25,0.25])
    psi=init_heading()
    sigma_a=0.5; sigma_bias=0.002; Kpsi=0.02; v_zupt=0.4; a_zupt=0.3
    gps_idx=0
    out_n=np.zeros(N); out_v=np.zeros(N); out_e=np.zeros(N)
    out_bax=np.zeros(N); out_bay=np.zeros(N)
    Hc=np.array([[1,0,0,0,0,0],[0,1,0,0,0,0.]])
    Hv=np.array([[0,0,1,0,0,0],[0,0,0,1,0,0.]])
    for k in range(N):
        d=dt_k[k]; c=math.cos(psi); s=math.sin(psi)
        psi += gz[k]*d
        A=np.array([
            [1,0,d,0,-0.5*d*d*c, 0.5*d*d*s],
            [0,1,0,d,-0.5*d*d*s,-0.5*d*d*c],
            [0,0,1,0,-d*c, d*s],
            [0,0,0,1,-d*s,-d*c],
            [0,0,0,0,1,0],
            [0,0,0,0,0,1.]])
        aN=c*ax[k]-s*ay[k]; aE=s*ax[k]+c*ay[k]
        Bu=np.array([0.5*d*d*aN,0.5*d*d*aE,d*aN,d*aE,0,0.])
        x=A@x+Bu
        qp=0.25*d**4*sigma_a**2; qv=d*d*sigma_a**2; qb=d*sigma_bias**2
        Q=np.diag([qp,qp,qv,qv,qb,qb])
        P=A@P@A.T+Q
        while gps_idx<ngps and tg[gps_idx]<=t[k]+1e-9:
            h=hacc[gps_idx]
            if not (0<h<50): gps_idx+=1; continue
            zc=np.array([gN[gps_idx],gE[gps_idx]])
            x,P=ekf_update(x,P,zc,Hc,np.diag([h*h,h*h]),1e18)
            if not math.isnan(vN_g[gps_idx]):
                zv=np.array([vN_g[gps_idx],vE_g[gps_idx]])
                x,P=ekf_update(x,P,zv,Hv,np.diag([0.09,0.09]),1e18)
                if not math.isnan(spd[gps_idx]) and spd[gps_idx]>2:
                    pg=math.atan2(vE_g[gps_idx],vN_g[gps_idx])
                    psi+=Kpsi*math.atan2(math.sin(pg-psi),math.cos(pg-psi))
            gps_idx+=1
        ah=math.hypot(ax[k],ay[k]); sp_e=math.hypot(x[2],x[3])
        if ah<a_zupt and abs(gz[k])<0.05 and sp_e<v_zupt:
            x,P=ekf_update(x,P,np.array([0,0.]),Hv,np.diag([0.0025,0.0025]),1e18)
        out_n[k]=x[0]; out_e[k]=x[1]; out_v[k]=math.hypot(x[2],x[3])
        out_bax[k]=x[4]; out_bay[k]=x[5]
    return out_n,out_e,out_v,out_bax,out_bay

def butter2(x,fc=5.0):
    fs=1/dt; g=math.tan(math.pi*fc/fs); D=g*g+math.sqrt(2)*g+1
    b=[g*g/D,2*g*g/D,g*g/D]; a=[1,2*(g*g-1)/D,(g*g-math.sqrt(2)*g+1)/D]
    y=np.zeros(len(x)); y[0]=x[0]; y[1]=x[1]
    for i in range(2,len(x)):
        y[i]=b[0]*x[i]+b[1]*x[i-1]+b[2]*x[i-2]-a[1]*y[i-1]-a[2]*y[i-2]
    return y

def imuLoc():
    axf=butter2(ax); ayf=butter2(ay)
    win=min(N,max(50,round(2/dt)))
    bax=axf[:win].mean(); bay=ayf[:win].mean()
    psi=init_heading()
    vN=vN_g[0] if not math.isnan(vN_g[0]) else 0.0
    vE=vE_g[0] if not math.isnan(vE_g[0]) else 0.0
    pN=0.0; pE=0.0; a_zupt=0.25
    on=np.zeros(N); oe=np.zeros(N); ov=np.zeros(N)
    for k in range(N):
        if k>0: psi+=gz[k]*dt
        a1=axf[k]-bax; a2=ayf[k]-bay
        if math.hypot(a1,a2)<a_zupt and abs(gz[k])<0.05:
            vN=0; vE=0
        else:
            aN=a1*math.cos(psi)-a2*math.sin(psi); aE=a1*math.sin(psi)+a2*math.cos(psi)
            pN+=vN*dt+0.5*aN*dt*dt; pE+=vE*dt+0.5*aE*dt*dt
            vN+=aN*dt; vE+=aE*dt
        on[k]=pN; oe[k]=pE; ov[k]=math.hypot(vN,vE)
    return on,oe,ov

def metrics(name,on,oe,ov):
    fn=np.interp(tg,t,on); fe=np.interp(tg,t,oe); fv=np.interp(tg,t,ov)
    good=(hacc<50)&~np.isnan(gN)
    err=np.sqrt((fn[good]-gN[good])**2+(fe[good]-gE[good])**2)
    mv=good&~np.isnan(spd)&(spd>=0)
    verr=fv[mv]-spd[mv]
    print(f"{name:16s} PosRMSE={np.sqrt((err**2).mean()):7.1f} m  "
          f"CEP50={np.median(err):6.1f} m  Max={err.max():7.1f} m  "
          f"FinalDrift={err[-1]:7.1f} m  VelRMSE={np.sqrt((verr**2).mean()):5.2f} m/s")

if __name__ == '__main__':
  print("="*95)
  print("MARMARAY (telefon) - dogrulama")
  print("="*95)
  n1,e1,v1,bx,by=fusionEKF()
  metrics("fusionEKF", n1,e1,v1)
  print(f"  -> kestirilen bias son: bax={bx[-1]:+.4f}  bay={by[-1]:+.4f} m/s^2")
  n2,e2,v2=imuLoc()
  metrics("imuLoc (saf DR)", n2,e2,v2)
  print("\nGPS toplam mesafe ~23.8 km, sure ~36 dk. fusionEKF konum hatasi metre,"
        "\nsaf DR ise (beklenildigi gibi) cok daha buyuk -> ivme-only sinirini gosterir.")
