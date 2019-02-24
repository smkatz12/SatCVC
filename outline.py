import numpy as np

main():
    figureOutEnv()
    Lloyd()

figureOutEnv():
    computePhi()

Lloyd(pInit, phi):
    computeVi()
    computeCVi()
    update_pi()
    if pi==CVi:
        end

arcDistance(p,q,R):
#Comment more

computePhi():
    do this everywhere
    normal = computeNormal(theta,phi)
    phi = abs(normal.dot(linetocenter))

computeNormal(theta,phi):

computeVi():
    use arcDistance(p,q,R)
    should return arrays of [theta, psi]

computeCVi():
    call computePhi()

update_pi(): 
    pidot = K*(Cvi-pi)
    pi_new = pi + pidot*dt 

