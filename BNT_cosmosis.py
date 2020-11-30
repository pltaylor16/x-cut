import numpy as np 
from cosmosis.datablock import names, option_section
from BNT import BNT
from math import pi as pi


def setup(options):

    BNT_path = options.get(option_section, 'BNT_path')
    run_mode = options.get(option_section, 'run_mode')
    BNT_matrix = np.loadtxt(BNT_path)
    

    return BNT_matrix, run_mode



def execute(block, config):

    BNT_matrix, run_mode = config

    nbins = BNT_matrix.shape[0]
 
    #k-cut in harmonic space
    if run_mode == 'harmonic':
        #load and format data from the bloc
        ell = block['shear_cl', 'ell'] 
        ell_shape = ell.shape[0]
        cl_data = np.zeros((nbins, nbins, ell_shape))
        for i in range(nbins):
            for j in range(nbins):
                if i >= j:
                    cl_data[i,j,:] = block['shear_cl', 'bin_%s_%s' %(i+1, j+1)] 
                    cl_data[j,i,:] = cl_data[i,j,:]


        #make the BNT transformation 
        cl_data = np.zeros((nbins, nbins, ell_shape))
        cl_data_new = np.zeros((nbins, nbins, ell_shape))
        for l in range(ell_shape):
            cl_data_new[:,:,l] = np.dot(np.dot(BNT_matrix, cl_data[:,:,l]), BNT_matrix.T)


        #save back to the block
        for i in range(nbins):
            for j in range(nbins):
                if i >= j:
                    block['shear_cl', 'bin_%s_%s' %(i+1, j+1)]  = cl_data_new[i,j,:] 


    #x-cut in configuration space
    elif run_mode == 'configuration':
        
        #load and format data from the bloc
        try:
            fmat = 'new'
            theta = block['shear_xi_plus', 'theta']
        except:
            fmat = 'old'
            theta = block['shear_xi', 'theta']

        theta_shape = theta.shape[0]
        xi_p_data = np.zeros((nbins, nbins, theta_shape))
        xi_m_data = np.zeros((nbins, nbins, theta_shape))
        xi_p_data_new = np.zeros((nbins, nbins, theta_shape))
        xi_m_data_new = np.zeros((nbins, nbins, theta_shape))

        for i in range(nbins):
            for j in range(nbins):
                if i >= j:
                    if fmat == 'new':
                        xi_p_data[i,j,:] = block['shear_xi_plus', 'bin_%s_%s' %(i+1, j+1)] 
                        xi_p_data[j,i,:] = xi_p_data[i,j,:]
                        xi_m_data[i,j,:] = block['shear_xi_minus', 'bin_%s_%s' %(i+1, j+1)] 
                        xi_m_data[j,i,:] = xi_m_data[i,j,:]
                    elif fmat == 'old':
                        xi_p_data[i,j,:] = block['shear_xi', 'xiplus_%s_%s' %(i+1, j+1)] 
                        xi_p_data[j,i,:] = xi_p_data[i,j,:]
                        xi_m_data[i,j,:] = block['shear_xi', 'ximinus_%s_%s' %(i+1, j+1)] 
                        xi_m_data[j,i,:] = xi_m_data[i,j,:]


        #make the BNT transformation and cut
        for t in range(theta_shape):
            xi_p_data_new[:,:,t] = np.dot(np.dot(BNT_matrix, xi_p_data[:,:,t]), BNT_matrix.T)
            xi_m_data_new[:,:,t] = np.dot(np.dot(BNT_matrix, xi_m_data[:,:,t]), BNT_matrix.T)


        #save back to the block
        for i in range(nbins):
            for j in range(nbins):
                if i >= j:
                    if fmat == 'new':
                        block['shear_xi_plus', 'bin_%s_%s' %(i+1, j+1)] = xi_p_data_new[i,j,:] 
                        block['shear_xi_minus', 'bin_%s_%s' %(i+1, j+1)] = xi_m_data_new[i,j,:] 
                    elif fmat == 'old':
                        block['shear_xi', 'xiplus_%s_%s' %(i+1, j+1)] = xi_p_data_new[i,j,:] 
                        block['shear_xi', 'ximinus_%s_%s' %(i+1, j+1)] = xi_m_data_new[i,j,:] 



        

    else:
        print ('run-mode must be either harmonic or configuration')

    return 0.



def cleanup(config):
    return 0.