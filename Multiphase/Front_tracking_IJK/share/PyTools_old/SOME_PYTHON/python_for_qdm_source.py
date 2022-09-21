# voir IJK_MEMO2 dans phone pour savoir ce que c'est ce fichier
import matplotlib.pyplot as plt
import numpy as np
import sys
import DNSTools3 as dtool

# emploi python python_for_qdm_source.py fichier_out nb_iteration_a_ignorer
# ATTENTION : pour utiliser ce script il faut etre au niveau des dossiers RUN0X

# DONNEES UTILISATEUR
chemin_acc=sys.argv[1]
start=int(sys.argv[2])
schema_temps=sys.argv[3]

# CHARGEMENT DONNEES
# x[0] -> qdm, x[1] -> uv, x[2] -> ul
x,t,i = np.loadtxt(chemin_acc)[:,-3:], np.loadtxt(chemin_acc)[:,1], np.loadtxt(chemin_acc)[:,0]

# Si simu RK3, on a les valeurs aux trois sous pas de temps (dont on se moque)
# -> On va donc supprimer les valeurs du premier et deuxieme sous pas de temps
# ATTENTION SI L'ECRITURE DES DONNES N'EST PLUS FAIT A CHAQUE RK3 PDT, 
# IL NE FAUT PLUS RENTRER DANS LA CONDITION DE SES MORTS CI-DESSOUS.
if schema_temps=="RK3":
    # Recupere le numero de run. Les noms des RUN sont : INIT, RUN00, RUN01, ...
    if "INIT" in chemin_acc:
        nrun=0
    else:
        nrun=int(chemin_acc.split("RUN")[-1][:2])+2
    # On a les trois sous pas de temps pour toutes les iterations sauf l'it 0
    i=np.insert(i,0,0)
    i=np.insert(i,0,0)
    # Chaque ligne correspond aux it d'un RUN. On reshape pour vectoriser le selecteur
    if nrun!=0:
        i=i.reshape(nrun,int((len(i))/nrun))
        # Pour chaque RUN enregistre, on ne garde qu'une occurence de l'iteration
        uniq_i,indices = np.unique(i,return_index=True,axis=1)
    else:
        # Pour chaque RUN enregistre, on ne garde qu'une occurence de l'iteration
        uniq_i,indices = np.unique(i,return_index=True)

    # Remise a plat de la liste des iterations epurees
    i=uniq_i.flatten()
    # La liste indices ne permet de recuperer que pour le premier RUN -le INIT- . 
    # -> Il faut concatener pour obtenir la liste COMPLETE des indices a travers les RUN.
    indices_tot=indices
    for r in range(nrun-1):
        indices_tot=np.concatenate((indices_tot,indices_tot[-len(indices):]+indices[-1]+3))
    # Selection des donnees du dernier sous pas de temps.
    x,t=x[indices_tot,:],t[indices_tot]

jdd = chemin_acc.split('/')[-1].split('_acceleration')[0]
jdd = chemin_acc.split('RUN')[0]+'RUN'+chemin_acc.split('RUN')[-1][:2]+'/'+jdd+'.data'
# ATTENTION SI L'ECRITURE DES DONNES N'EST PLUS FAIT A CHAQUE RK3 PDT, 
# IL NE FAUT PLUS RENTRER DANS LA CONDITION DE SES MORTS CI-DESSUS.

# Selection partielle
xs,ts,ist = x[start:,:],t[start:],i[start:]
#####################################################################################
# GRAPHIQUE I : qdm - uv 

fig, (ax1,ax2) = plt.subplots(2,1)
# ins1, ins2  = ax1.inset_axes([0.3, 0.2, 0.6, 0.6]), ax2.inset_axes([0.3, 0.2, 0.6, 0.6])

ax1.plot(t,x[:,1],'k',label='u_v');
ax1.plot(t[np.argwhere(i==0)][:,0],x[np.argwhere(i==0),1],'oy');
# ins1.plot(ts,xs[:,1],'k')
# ins1.plot(ts[np.argwhere(ist==0)][:,0],xs[np.argwhere(ist==0),1],'oy');
# ins1.tick_params(axis="y",direction="in")#,pad=-0)
# ins1.tick_params(axis="x",direction="in")#,pad=-15)
ax1.legend();

ax2.plot(t,x[:,0],'r',label='qdm');
ax2.plot(t[np.argwhere(i==0)][:,0],x[np.argwhere(i==0),0],'og');
# ins2.plot(ts,xs[:,0],'r')
# ins2.plot(ts[np.argwhere(ist==0)][:,0],xs[np.argwhere(ist==0),0],'og');
# ins2.tick_params(axis="y",direction="in")#,pad=-0)
# ins2.tick_params(axis="x",direction="in")#,pad=-15)
ax2.set_xlabel("temps (s)")
ax2.legend();

for ax in (ax1,ax2):
    ax.tick_params(labelsize=22)
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
            label.set_fontsize(22)

fig.tight_layout();
fig.savefig("qdm_uv.png")

del(fig,ax1,ax2)#,ins1,ins2)
########################################################################################
# GRAPHIQUE II : uv - ul

fig1, (ax11,ax21) = plt.subplots(2,1)
ins11, ins21  = ax11.inset_axes([0.3, 0.2, 0.6, 0.6]), ax21.inset_axes([0.3, 0.2, 0.6, 0.6])

ax11.plot(t,x[:,1],'k',label='u_v');
ax11.plot(t[np.argwhere(i==0)][:,0],x[np.argwhere(i==0),1],'oy');
ins11.plot(ts,xs[:,1],'k')
ins11.plot(ts[np.argwhere(ist==0)][:,0],xs[np.argwhere(ist==0),1],'oy');
ins11.tick_params(axis="y",direction="in")#,pad=-0)
ins11.tick_params(axis="x",direction="in")#,pad=-15)
ax11.legend();

ax21.plot(t,x[:,2],'b',label='u_l');
ax21.plot(t[np.argwhere(i==0)][:,0],x[np.argwhere(i==0),2],'om');
ins21.plot(ts,xs[:,2],'b')
ins21.plot(ts[np.argwhere(ist==0)][:,0],xs[np.argwhere(ist==0),2],'om');
ins21.tick_params(axis="y",direction="in")#,pad=-0)
ins21.tick_params(axis="x",direction="in")#,pad=-15)
ax21.set_xlabel("temps (s)")
ax21.legend();

for ax in (ax11,ax21):
    ax.tick_params(labelsize=22)
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
            label.set_fontsize(22)

fig1.tight_layout();
fig1.savefig("ul_uv.png")

########################################################################################
# GRAPHIQUE III : ur=uv-ul - qdm-ul
d = x[:,1] - x[:,2]
ds = xs[:,1] - xs[:,2]
e = x[:,0] - x[:,2]
es = xs[:,0] - xs[:,2]


fig2, (ax12,ax22) = plt.subplots(2,1)
ins12, ins22  = ax12.inset_axes([0.3, 0.2, 0.6, 0.6]), ax22.inset_axes([0.3, 0.2, 0.6, 0.6])

ax12.plot(t,d,'g',label='u_r');
ax12.plot(t[np.argwhere(i==0)][:,0],d[np.argwhere(i==0)][:,0],'ob');
ins12.plot(ts,ds,'g')
ins12.plot(ts[np.argwhere(ist==0)][:,0],ds[np.argwhere(ist==0)][:,0],'ob');
ins12.tick_params(axis="y",direction="in")#,pad=-0)
ins12.tick_params(axis="x",direction="in")#,pad=-15)
ax12.legend();

ax22.plot(t,e,'y',label='qdm-u_l');
ax22.plot(t[np.argwhere(i==0)][:,0],e[np.argwhere(i==0)][:,0],'ok');
ins22.plot(ts,es,'y')
ins22.plot(ts[np.argwhere(ist==0)][:,0],es[np.argwhere(ist==0)][:,0],'ok');
ins22.tick_params(axis="y",direction="in")#,pad=-0)
ins22.tick_params(axis="x",direction="in")#,pad=-15)
ax22.set_xlabel("temps (s)")
ax22.legend();

for ax in (ax12,ax22):
    ax.tick_params(labelsize=22)
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
            label.set_fontsize(22)

fig2.tight_layout();
fig2.savefig("ur_ul-qdm.png")

########################################################################################
# GRAPHIQUE IV : qdm adimentionne : qdm/(terme_force_init) = qdm/(rho_m g) et qdm/(u_v')
terme_force_init = float(dtool.getParam(jdd, "terme_force_init"))
if terme_force_init != 0:
    qdm_adim0 = x[:,0] / terme_force_init 
    qdm_adim0s = xs[:,0] / terme_force_init 
# qdm_adim0 = x[:,0] / dtool.getParam("terme_force_init",jdd) 


fig2, (ax12,ax22) = plt.subplots(2,1)
# ins12, ins22  = ax12.inset_axes([0.3, 0.2, 0.6, 0.6]), ax22.inset_axes([0.3, 0.2, 0.6, 0.6])

ax12.plot(t,qdm_adim0,'g',label=r'qdm/$\overline{\rho} g$');
ax12.plot(t[np.argwhere(i==0)][:,0],qdm_adim0[np.argwhere(i==0)][:,0],'ob');
# ins12.plot(ts,qdm_adim0s,'g')
# ins12.plot(ts[np.argwhere(ist==0)][:,0],qdm_adim0s[np.argwhere(ist==0)][:,0],'ob');
# ins12.tick_params(axis="y",direction="in")#,pad=-0)
# ins12.tick_params(axis="x",direction="in")#,pad=-15)
ax12.legend();

# ax22.plot(t,e,'y',label='qdm-u_l');
# ax22.plot(t[np.argwhere(i==0)][:,0],e[np.argwhere(i==0)][:,0],'ok');
# ins22.plot(ts,es,'y')
# ins22.plot(ts[np.argwhere(ist==0)][:,0],es[np.argwhere(ist==0)][:,0],'ok');
# ins22.tick_params(axis="y",direction="in")#,pad=-0)
# ins22.tick_params(axis="x",direction="in")#,pad=-15)
# ax22.set_xlabel("temps (s)")
# ax22.legend();

for ax in (ax12,ax22):
    ax.tick_params(labelsize=22)
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
            label.set_fontsize(22)

fig2.tight_layout();
fig2.savefig("qdm_div_par_terme_force_init.png")


# MONTRE TOUS
plt.show()

