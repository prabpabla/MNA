# -*- coding: utf-8 -*-
import numpy as np
import scipy
import scipy.linalg
import matplotlib.pyplot as plt

class MatrixConnection:
    
    # Single matrix connection data
    def __init__(self, node1=-1, node2=-1, sym="Z", value=-1.0):
        self.node1 = node1
        self.node2 = node2
        self.sym = sym
        self.value = value
    
    def get_node1(self):
        return self.node1
    
    def get_node2(self):
        return self.node2
    
    def get_sym(self):
        return self.sym
          
    def get_value(self):
        return self.value
        
    def is_res(self):
        return self.sym == "G"
        
    def is_cap(self):
        return self.sym == "C"
    
    def is_ind(self):
        return self.sym == "L"
    
    def is_cur(self):
        return self.sym == "S"
    
    def is_vol(self):
        return self.sym == "V"


    
if __name__ == '__main__':
    file_name = "netlist2.txt"
    
    ind_sym = "L"
    vol_sym = "V"
    
    matrix_size = 0
    matrix_size_before_ind = 0
    ind_count = 0
    vol_count = 0
    matrix_connections = []
    h = 0.01e-9 #time step
    time = 0.0
    itteration = 3000
    
    with open(file_name) as f:
        for line in f:  
            numbers_str = line.split() #split line data
            
            #values from file
            node1 = int(numbers_str[0])
            node2 = int(numbers_str[1])
            sym = numbers_str[2]
            value = float(numbers_str[3])
            
            connection = MatrixConnection(node1,node2,sym,value)
            matrix_connections.append(connection)       
            
            #get largest node
            if(node1>matrix_size):
                matrix_size = node1
            elif(node2>matrix_size):
                matrix_size = node2
            
            #get indcutor count
            if(sym == ind_sym):
                ind_count = ind_count + 1
            
            #get voltage count
            if(sym == vol_sym):
                vol_count = vol_count + 1
    
    # matrix_size_before_ind = matrix_size
    # matrix_size = int(matrix_size + ind_count + vol_count)
    G_matrix = np.zeros((matrix_size,matrix_size))
    C_matrix = np.zeros((matrix_size,matrix_size))
    B_matrix = np.zeros((matrix_size,1))

    #B_matrix  = np.resize(B_matrix,(matrix_size+1,1))  
    #B_matrix[matrix_size][0] = 0
    
    print("Matrix size will be ",matrix_size)
    
    # update matrix for all resistors in the connection list
    for i in range(0,len(matrix_connections)):
        
        connection = matrix_connections[i]
        node1 = connection.get_node1()
        node2 = connection.get_node2()
        value = connection.get_value()
        
        # C and G Matrix operations
        if(connection.is_res()): #if symbol is resistor 
        
            if(not (node1 == 0)):
                # node 1 is not zero
                G_matrix[node1-1][node1-1] += 1.0/value                
                
            if(not (node2 == 0)):
                # node 2 is not zero
                G_matrix[node2-1][node2-1] += 1.0/value
                
            if(not (node1 == 0) and not (node2 == 0)):
                # both nodes are not zero
                G_matrix[node1-1][node2-1] += -1.0/value
                G_matrix[node2-1][node1-1] += -1.0/value
  
        elif(connection.is_cap()):  #if symbol is capacitor
            
            if(not (node1 == 0)):
                # node 1 is not zero
                C_matrix[node1-1][node1-1] += value                
                
            if(not (node2 == 0)):
                # node 2 is not zero
                C_matrix[node2-1][node2-1] += value
                
            if(not (node1 == 0) and not (node2 == 0)):
                # both nodes are not zero
                C_matrix[node1-1][node2-1] += -value
                C_matrix[node2-1][node1-1] += -value  
                
        elif(matrix_connections[i].is_ind()):  #if symbol is inductor 

            result = np.zeros((matrix_size+1,matrix_size+1))
            result[:G_matrix.shape[0],:G_matrix.shape[1]] = G_matrix
            G_matrix = result
            
            result = np.zeros((matrix_size+1,matrix_size+1))
            result[:C_matrix.shape[0],:C_matrix.shape[1]] = C_matrix
            C_matrix = result

            B_matrix  = np.resize(B_matrix,(matrix_size+1,1))  
            B_matrix[matrix_size][0] = 0

            matrix_size = matrix_size + 1
            
            if(not (node1 == 0)):
                # node 1 is not zero
                G_matrix[node1-1][matrix_size-1] += 1.0              
                G_matrix[matrix_size-1][node1-1] += 1.0        
                
            if(not (node2 == 0)):
                # node 2 is not zero
                G_matrix[node2-1][matrix_size-1] += -1.0               
                G_matrix[matrix_size-1][node2-1] += -1.0                   
                
            if(not (node1 == 0) and not (node2 == 0)):
                # both nodes are not zero
                C_matrix[matrix_size-1][matrix_size-1] += -value
        
        elif(matrix_connections[i].is_vol()):  #if symbol is voltage
            result = np.zeros((matrix_size+1,matrix_size+1))
            result[:G_matrix.shape[0],:G_matrix.shape[1]] = G_matrix
            G_matrix = result
            
            result = np.zeros((matrix_size+1,matrix_size+1))
            result[:C_matrix.shape[0],:C_matrix.shape[1]] = C_matrix
            C_matrix = result
            
            B_matrix  = np.resize(B_matrix,(matrix_size+1,1))  
            B_matrix[matrix_size][0] = 0
            
            matrix_size = matrix_size + 1
            
            B_matrix[matrix_size-1][0] = value
             
            if(not (node1 == 0)):
                # node 1 is not zero
                G_matrix[node1-1][matrix_size-1] += 1.0              
                G_matrix[matrix_size-1][node1-1] += 1.0        
                
            if(not (node2 == 0)):
                # node 2 is not zero
                G_matrix[node2-1][matrix_size-1] += -1.0               
                G_matrix[matrix_size-1][node2-1] += -1.0               
                
                          
                       
        
    print("Matrix size will be ",matrix_size)    
    print("G Matrix:\n",G_matrix) #print G matrix
    print("-----------------------------------")
    print("C Matrix:\n",C_matrix) #print C matrix
    print("-----------------------------------")
    print("B Matrix:\n",B_matrix) #print B matrix
    print("-----------------------------------")
    
    
    ######################################################################
    
    
    C_div_H_matrix = np.zeros((matrix_size,matrix_size))
    #G_plus_C_div_H_matrix = np.zeros((matrix_size,matrix_size))
    X_matrix = np.zeros((matrix_size,1))
    X_matrix_storage = []
    
    
    C_div_H_matrix = np.divide(C_matrix,h)
    time_input = []
    
    
    time = 0.0
    value = 0.0
    time_plot = []
    # Matrixies that are time dependent 
    for i in range(0,itteration):
        
        #temp_B_matrix = np.copy(B_matrix)
        temp_B_matrix = np.zeros((matrix_size,1))
        for j in range(0,matrix_size):
            temp_B_matrix[j][0] = B_matrix[j][0]   
        
        
        if time<=2e-9:
            value=0.0
        elif time>2e-9 and time<=3e-9:
            value=(5/1e-9) * time - 10
        elif time>3e-9 and  time <=13e-9:
            value=5.0
        elif time > 13e-9 and time <= 14e-9:
            value=((-5/1e-9)) * time + 70
        else: 
            value=0.0
           
        # this for loop only gets run if there is a voltage source
        for j in range(0,len(matrix_connections)):
        
            if(matrix_connections[j].is_cur()):  #if symbol is current source 
            
                if(not (node1 == 0)):
                    temp_B_matrix[node1 - 1] = -value
                    
                if(not (node2 == 0)):
                    temp_B_matrix[node2 - 1] = value
        
        
        temp_B_matrix = np.dot(temp_B_matrix,value)
        
        time_input.append(temp_B_matrix)
        time_plot.append(value)     #value) 

        time = time + h
    
    plt.plot(time_plot)
    plt.show()
     
    
    ''' 
    # My solution to X matrix
    G_plus_C_div_H_matrix = np.add(C_div_H_matrix,G_matrix) # A
    X_prev_times_C_div_H_plus_B_matrix = np.add(np.dot(C_div_H_matrix,X_prev_matrix),B_matrix) # BN
    G_plus_C_div_H_matrix_transpose = np.transpose(G_plus_C_div_H_matrix)
    X_matrix = np.dot(G_plus_C_div_H_matrix_transpose,X_prev_times_C_div_H_plus_B_matrix)
    ''' 
    time = 0.0
    A = np.add(C_div_H_matrix,G_matrix) # A
    P,L,U = scipy.linalg.lu(A)
    #print(L)
    
    for i in range(0,itteration-1):
        #B_matrix[7][0] = time_plot[i+1]
        #G_plus_C_div_H_matrix = np.add(C_div_H_matrix,G_matrix) # A
        #X_prev_times_C_div_H_plus_B_matrix = np.add(np.dot(C_div_H_matrix,X_matrix),B_matrix) # BN
        #G_plus_C_div_H_matrix_transpose = np.transpose(G_plus_C_div_H_matrix)
        #X_matrix = np.dot(G_plus_C_div_H_matrix_transpose,X_prev_times_C_div_H_plus_B_matrix)
        
        BN = np.add(np.dot(C_div_H_matrix,X_matrix), time_input[i+1]) # left side

        Z = np.linalg.solve(L,(np.dot(P,BN)))
        X_matrix = np.linalg.solve(U,Z)
        
        #BOTH ARE THE SAME, BOTH WORK
        #Z = np.dot(np.linalg.inv(L),(np.dot(P,BN)))
        #X_matrix = np.dot(np.linalg.inv(U),Z)
        
        X_matrix_storage.append(X_matrix[6][0])
        
        time = time + h
        
        #if(i>=205 and i<=207):
            #print(X_matrix)
            
        
        
    plt.plot(X_matrix_storage)
    plt.show()
    
    #print("B Matrix:\n",B_matrix) #print B matrix
    #print("-----------------------------------")
    '''
    print("C div H Matrix:\n",C_div_H_matrix) #print B matrix
    print("-----------------------------------")
    print("G plus C div H Matrix:\n",A) 
    print("-----------------------------------")
    print("X prev times C div H plus B Matrix:\n",BN) 
    print("-----------------------------------")
    print("X Matrix:\n",X_matrix) 
    '''
    '''
    1 0 G 1 0.0001
    1 2 G 2 0.0005
    3 0 G 3 0.0002
    2 0 C 1 0.0000000001
    2 3 L 1 0.0000000001
    0 1 S 1 0.0001
    '''