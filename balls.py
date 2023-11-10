#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

############  2 Balls functions  ################
# RO is a point
# RQ and DQO are arrays of same length
# Most of those functions are readable even with RO RQ and DQO of same shape

# two balls arent apprt?
def isnt_appart(R1,R2,D):
    return (D<(R1+R2))


# one ball is contained in the other?
def Is_inside(RO,RQ,DQO):
    return RO+DQO<RQ, RQ+DQO<RO

#Volume of a ball    
def vol_1(R):
    return 4*np.pi*R**3/3

# Tthe volume of a lense (non-trivial intersection of two balls) is given by
def co_volume_2(R1,R2,D):    
    return np.pi*((R1+R2-D)**2)*(D**2-3*(R1-R2)**2+2*D*(R1+R2))/12/D


    

# combiine all cases above    
def co_vol_2(RO,RQ,DQO):
    # keep those who are not appart
    QO_isnt_apt_inds=isnt_appart(RO,RQ,DQO);
    RQt=RQ[QO_isnt_apt_inds]
    DQOt=DQO[QO_isnt_apt_inds]
    
    vol_2_temp=np.zeros(len(RQt))

    
    # find simple cases and fill array
    O_in_Q_inds, Q_in_O_inds= Is_inside(RO,RQt,DQOt)
    vol_2_temp[O_in_Q_inds]=np.full(np.sum(O_in_Q_inds),vol_1(RO))
    vol_2_temp[Q_in_O_inds]=vol_1(RQt[Q_in_O_inds])
    
    # find complicated cases
    QO_not_in=np.logical_not((O_in_Q_inds) + (Q_in_O_inds))
    vol_2_temp[QO_not_in]=co_volume_2(RO,RQt[QO_not_in],DQOt[QO_not_in])
    return vol_2_temp, QO_isnt_apt_inds, O_in_Q_inds, Q_in_O_inds, QO_not_in, RQt, DQOt


############  3 Balls functions  ################

######### Some auxilary variables
def eps_q_w(R1,R2,R3,D1,D2,D3):
    R1s=R1**2
    R2s=R2**2
    R3s=R3**2
    D1s=D1**2
    D2s=D2**2
    D3s=D3**2
    
    eps1=(R2s-R3s)/D1s
    eps2=(R3s-R1s)/D2s
    eps3=(R1s-R2s)/D3s
    
    q1=D1*(D2s+D3s-D1s+R2s+R3s-2*R1s+eps1*(D2s-D3s))
    q2=D2*(D3s+D1s-D2s+R3s+R1s-2*R2s+eps2*(D3s-D1s))
    q3=D3*(D1s+D2s-D3s+R1s+R2s-2*R3s+eps3*(D1s-D2s))
    
    w2=(R1s*D1s+R2s*D2s+R3s*D3s)*(D1s+D2s+D3s)-2*(R1s*D1s**2+R2s*D2s**2+R3s*D3s**2)+D1s*D2s*D3s*(eps1*eps2+eps2*eps3+eps3*eps1-1)
    
    return eps1, eps2, eps3, q1, q2, q3, w2


######### The non-trivial case - topology d. 
######### A combination of my calculation and K. D. Gibson and A. Scheraga, Tech. Rep. (1987)
def vol_3_topo_d(R1,R2,R3,D1,D2,D3,eps1, eps2, eps3, q1, q2, q3, w2):
    R1s=R1**2
    R2s=R2**2
    R3s=R3**2
    D1s=D1**2
    D2s=D2**2
    D3s=D3**2
    
    R1q=R1**3
    R2q=R2**3
    R3q=R3**3
    
    # Volume of tetrahedron = w/12
    w=np.sqrt(w2)
    
    # At each vertex of the tetrahedron we define the three angles a_i,b_i,c_i
    a1=np.arccos((D2s+D3s-D1s)/(2*D2*D3))
    a2=np.arccos((D3s+D1s-D2s)/(2*D3*D1))
    a3=np.arccos((D1s+D2s-D3s)/(2*D1*D2))
    
    b1=np.arccos((R1s+D2s-R3s)/(2*R1*D2))
    b2=np.arccos((R2s+D3s-R1s)/(2*R2*D3))
    b3=np.arccos((R3s+D1s-R2s)/(2*R3*D1))
    
    g1=np.arccos((D3s+R1s-R2s)/(2*D3*R1))
    g2=np.arccos((D1s+R2s-R3s)/(2*D1*R2))
    g3=np.arccos((D2s+R3s-R1s)/(2*D2*R3))
    
    # The angles are than used to calculate the solid angle using l'Huilier's theorem.
    s1=(a1+b1+g1)/2
    s2=(a2+b2+g2)/2
    s3=(a3+b3+g3)/2
    # l'Huilier's theorem:
    O1=4*np.arctan2(np.sqrt(np.tan(s1/2)*np.tan((s1-a1)/2)*np.tan((s1-b1)/2)*np.tan((s1-g1)/2)),1)
    O2=4*np.arctan2(np.sqrt(np.tan(s2/2)*np.tan((s2-a2)/2)*np.tan((s2-b2)/2)*np.tan((s2-g2)/2)),1)
    O3=4*np.arctan2(np.sqrt(np.tan(s3/2)*np.tan((s3-a3)/2)*np.tan((s3-b3)/2)*np.tan((s3-g3)/2)),1)
    # Those angles are used to calculate (twice) the volume of each sphere that is bounded by the relevant faces of the tetrahedron
    T1=2*R1q*O1/3
    T2=2*R2q*O2/3
    T3=2*R3q*O3/3

    # The angle between relevant faces of the tetrahedron
    atan1=np.arctan2(2*w,q1)
    atan2=np.arctan2(2*w,q2)
    atan3=np.arctan2(2*w,q3)
 	# Those angles are used to calculate the relevant part of the intersection of two balls. 
    U1=atan1*co_volume_2(R2,R3,D1)/np.pi
    U2=atan2*co_volume_2(R3,R1,D2)/np.pi
    U3=atan3*co_volume_2(R1,R2,D3)/np.pi
    
    
    #and we return the relevant volume
    return w/6-T1-T2-T3+U1+U2+U3

######### The non-trivial case - topologies a-c. 
######### My intrpretation of K. D. Gibson and A. Scheraga, Tech. Rep. (1987)
def vol_3_topo_ac(R1,R2,R3,D1,D2,D3):
    vol_3_temp=np.zeros(len(R1))
    
    R1s=R1**2
    R2s=R2**2
    R3s=R3**2
    D1s=D1**2
    D2s=D2**2
    D3s=D3**2
    
    # Define semiperimeters
    S=(D1+D2+D3)/2
    S1=(D1+R2+R3)/2
    S2=(D2+R3+R1)/2
    S3=(D3+R1+R2)/2
    
    # Areas of relvant trianles
    A=np.sqrt(S*(S-D1)*(S-D2)*(S-D3))
    A1=np.sqrt(S1*(S1-D1)*(S1-R2)*(S1-R3))
    A2=np.sqrt(S2*(S2-D2)*(S2-R3)*(S2-R1))
    A3=np.sqrt(S3*(S3-D3)*(S3-R1)*(S3-R2))
    
    # and distances form sphere sphere intersections
    L1p=(D2s-D3s+R2s-R3s)**2+16*(A+A1)**2-4*D1s*R1s
    L1m=L1p-16*4*A*A1
    L2p=(D3s-D1s+R3s-R1s)**2+16*(A+A2)**2-4*D2s*R2s
    L2m=L2p-16*4*A*A2
    L3p=(D1s-D2s+R1s-R2s)**2+16*(A+A3)**2-4*D3s*R3s
    L3m=L3p-16*4*A*A3
    
    # The conditions
    L1p_inds=L1p>0
    L1m_inds=L1m>0
    L2p_inds=L2p>0
    L2m_inds=L2m>0
    L3p_inds=L3p>0
    L3m_inds=L3m>0
    
    # case 1
    all_but_1_pos=np.logical_not(L1p_inds+L1m_inds)*L2p_inds*L2m_inds*L3p_inds*L3m_inds
    all_but_2_pos=np.logical_not(L2p_inds+L2m_inds)*L3p_inds*L3m_inds*L1p_inds*L1m_inds
    all_but_3_pos=np.logical_not(L3p_inds+L3m_inds)*L1p_inds*L1m_inds*L2p_inds*L2m_inds
    
    # case 2
    all_but_1_neg=L1p_inds*L1m_inds*np.logical_not(L2p_inds+L2m_inds)*np.logical_not(L3p_inds+L3m_inds)
    all_but_2_neg=L2p_inds*L2m_inds*np.logical_not(L3p_inds+L3m_inds)*np.logical_not(L1p_inds+L1m_inds)
    all_but_3_neg=L3p_inds*L3m_inds*np.logical_not(L1p_inds+L1m_inds)*np.logical_not(L2p_inds+L2m_inds)
    
    # case 1 fill
    vol_3_temp[all_but_1_pos]=co_volume_2(R2[all_but_1_pos],R3[all_but_1_pos],D1[all_but_1_pos])
    vol_3_temp[all_but_2_pos]=co_volume_2(R3[all_but_2_pos],R1[all_but_2_pos],D2[all_but_2_pos])
    vol_3_temp[all_but_3_pos]=co_volume_2(R1[all_but_3_pos],R2[all_but_3_pos],D3[all_but_3_pos])
    
    # case 2 fill
    Rit=R1[all_but_1_neg]
    Rjt=R2[all_but_1_neg]
    Rkt=R3[all_but_1_neg]
    Djt=D2[all_but_1_neg]
    Dkt=D3[all_but_1_neg]
    vol_3_temp[all_but_1_neg]=co_volume_2(Rit,Rjt,Dkt)+co_volume_2(Rkt,Rit,Djt)-vol_1(Rit)
    Rit=R2[all_but_2_neg]
    Rjt=R3[all_but_2_neg]
    Rkt=R1[all_but_2_neg]
    Djt=D3[all_but_2_neg]
    Dkt=D1[all_but_2_neg]
    vol_3_temp[all_but_2_neg]=co_volume_2(Rit,Rjt,Dkt)+co_volume_2(Rkt,Rit,Djt)-vol_1(Rit)
    Rit=R3[all_but_3_neg]
    Rjt=R1[all_but_3_neg]
    Rkt=R2[all_but_3_neg]
    Djt=D1[all_but_3_neg]
    Dkt=D2[all_but_3_neg]
    vol_3_temp[all_but_3_neg]=co_volume_2(Rit,Rjt,Dkt)+co_volume_2(Rkt,Rit,Djt)-vol_1(Rit)
    
    return vol_3_temp

######### The non-trivial cases combined
def vol_3_topo(R1,R2,R3,D1,D2,D3):
    eps1, eps2, eps3, q1, q2, q3, w2 =eps_q_w(R1,R2,R3,D1,D2,D3)
    vol_3_temp=np.zeros(len(R1))
    
    # Separate topology d from the rest
    case_d_inds=w2>=0
    not_d_inds=np.logical_not(case_d_inds)
    
    vol_3_temp[case_d_inds]=vol_3_topo_d(R1[case_d_inds],R2[case_d_inds],R3[case_d_inds],
                                         D1[case_d_inds],D2[case_d_inds],D3[case_d_inds],eps1[case_d_inds], eps2[case_d_inds], eps3[case_d_inds], q1[case_d_inds], q2[case_d_inds], q3[case_d_inds], w2[case_d_inds])
    
    vol_3_temp[not_d_inds]=vol_3_topo_ac(R1[not_d_inds],R2[not_d_inds],R3[not_d_inds],D1[not_d_inds],D2[not_d_inds],D3[not_d_inds])
    
    return vol_3_temp


######### A combination of all cases above the outputs both the 2-balls and 3-balls intersection
# RO 	- an integer: observer's radius.
# RQ 	- a list of length N: Quasars radii
# DQO 	- a list of length N: Quasar-Observer distances
# DQQ 	- array of size NxN : Quasar-Quasar distances
def co_vol_23(RO,RQ,DQO,DQQ):
    
    # Calculate 2-balls volume and relevant indecies
    vol_2_temp, QO_isnt_apt_inds, O_in_Q_inds, Q_in_O_inds, QO_not_in, RQt, DQOt = co_vol_2(RO,RQ,DQO)
    
    # Keep only relevant QQ pairs (QO aren't disjoint)
    DQQt=DQQ[QO_isnt_apt_inds][:,QO_isnt_apt_inds]
    
    # make a 2d arry out of vol_2
    vol_2_arr=np.full(DQQt.shape,vol_2_temp)
    
    # Find QQ pairs that aren't disjoint
    RQt_arr=np.full(DQQt.shape,RQt)
    QQ_isnt_apt_inds=isnt_appart(RQt_arr,np.transpose(RQt_arr),DQQt)
    

    
    # We will not sum over all i&j since we can sum over j>i and use our knowledge of N (c.f. eq. 8 in my note).
    # We therefore define an upper triangular boolean matrix (keeps only j>i)
    upper_tri_bol=np.triu(np.full(DQQt.shape,True),1)
    
    # initial empty volume arr
    vol_3_temp=np.zeros(DQQt.shape)

    
    ### Starting to fill array case by case
    
    # When RO ⊂ Qᵢ VO_{ij}=VO_j (and also i<-->j)
    vol_3_temp[O_in_Q_inds]=vol_2_arr[O_in_Q_inds]
    vol_3_temp[:,O_in_Q_inds]=np.transpose(vol_3_temp[O_in_Q_inds])
    # keep track of what still should be filled 
    O_in_Q_arr=np.full(DQQt.shape,O_in_Q_inds)
    not_O_in_Q_inds=np.logical_not(O_in_Q_arr+np.transpose(O_in_Q_arr))
    good_inds=not_O_in_Q_inds*QQ_isnt_apt_inds*upper_tri_bol
    
    # When Qⱼ ⊂ Qᵢ VO_{ij}=VO_j etc.
    # Find (relevant) QQ pairs that are inside one another
    temp_vol= vol_3_temp[good_inds]
    i_in_j_inds, j_in_i_inds= Is_inside(RQt_arr[good_inds],np.transpose(RQt_arr)[good_inds],DQQt[good_inds])
    temp_vol[i_in_j_inds]=vol_2_arr[good_inds][i_in_j_inds]
    temp_vol[j_in_i_inds]=np.transpose(vol_2_arr)[good_inds][j_in_i_inds]
    vol_3_temp[good_inds]=temp_vol
    # Update what's left to be filled
    temp_inds=good_inds[good_inds]
    temp_inds[j_in_i_inds+i_in_j_inds]=np.full(np.sum(j_in_i_inds+i_in_j_inds),False)
    good_inds[good_inds]=temp_inds
    
    # When [(Qᵢ⊂ RO)∩¬(RO⊂ Qⱼ)]∪[(Qⱼ⊂ RO)∩¬(RO⊂ Qᵢ)], we must calculate V_{ij}, the intersection of Qᵢ and Qⱼ
    Qi_in_O_not_in_Qj=np.outer(Q_in_O_inds,np.logical_not(O_in_Q_inds))
    Vij_bol_arr=(Qi_in_O_not_in_Qj+np.transpose(Qi_in_O_not_in_Qj))
    # and fill the relevant entries of the array
    vol_3_temp[Vij_bol_arr*good_inds]=co_volume_2(RQt_arr[Vij_bol_arr*good_inds],np.transpose(RQt_arr)[Vij_bol_arr*good_inds],DQQt[Vij_bol_arr*good_inds])
    # Update what's left to be filled
    good_inds[Vij_bol_arr]=False*good_inds[Vij_bol_arr]
    
    ## And to the complicated topologies
    long_R1=RQt_arr[good_inds]
    long_R2=np.transpose(RQt_arr)[good_inds]
    long_RO=0*long_R2+RO
    
    DQOt_arr=np.full(DQQt.shape,DQOt)
    long_DQ1=DQOt_arr[good_inds]
    long_DQ2=np.transpose(DQOt_arr)[good_inds]
    
    long_DQQ=DQQt[good_inds]
    
    vol_3_temp[good_inds]=vol_3_topo(long_RO,long_R1,long_R2,long_DQQ,long_DQ2,long_DQ1)
    
    
    return vol_3_temp[upper_tri_bol], vol_2_temp


