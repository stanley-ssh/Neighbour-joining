import sys


#To read a distance matrix file
#returns a list of taxa (list) and a distance matrix (dictionary, where the key is a pair of taxa and the value is the distance)
def readDistMatrixFile(filename):
    file1 = open(filename,"r")
    taxa_name = file1.readline().split()
    lines  = [] 
    # we read the ramining lines into a 2D list  (lines)
    for line in file1:
        lines .append( line.split())
        
    file1.close()

    # to assign we loop throgh the the 2D and store the distances in the dictionary 
    taxa_len = len(taxa_name)
    distance_dict = {} 
    

    for i  in range (taxa_len):
        # loop through 
        j=1 
        while ( j <= taxa_len):
            key = (taxa_name[i],taxa_name[j-1])
            distance_dict[key] = lines[i][j]

            j= j+1
    

    # Done 
    
    return (taxa_name,distance_dict)


    

#TO COMPLETE



#To produce the U table: union of distances for each taxon in the supplied list of taxa (also takes the distance matrix as a parameter)
#Returns the U table (a dictionary, where the key is the taxon and the value is the union of distances to all other taxa)
def produceUTable(listOfTaxa, distMatrix):
    #print (distMatrix)
    U_dict = {} 
    for taxa in listOfTaxa:
        sum_val = 0 
        for key in distMatrix:
            
            curr_taxa = key[0]

            if taxa == curr_taxa:

                sum_val = sum_val +  float(distMatrix[key]) 
        
        U_dict[taxa] = sum_val 

    return U_dict

     
#TO COMPLETE


#To produce a delta matrix from a distance matrix given as a parameter (also takes a list of taxa and the uTable as parameters).
#Returns the delta matrix (as a dictionary, where the key is a pair of taxa and the value is the delta value)
def produceDeltaMatrix(listOfTaxa, distMatrix, uTable):
    # formular : delta(i,j) = (N-2) Dij - Ui - Uj 
    N = len(listOfTaxa)

    delta_mat = {} 

    #print(listOfTaxa)
    i =0 
    j=0 
        
    for i in range(N):
        j = i + 1 
        while (j < N):
            value = listOfTaxa[i]
            val2 = listOfTaxa[j]
            key = (listOfTaxa[i],listOfTaxa[j])
            # print (" the key is {}".format(key))
            
            value = (N-2 ) * float(distMatrix[key]) - uTable[value] - uTable[val2]

            delta_mat[key] = value
            j+=1
    
    return delta_mat


     
#TO COMPLETE


#Recursive method that does all the steps of NJ, updates the distance matrix and calls itself on the updated distance matrix and list of taxa
def neighborJoiningRec(listOfTaxa, distMatrix):

    N= len(listOfTaxa)
    uTable   = produceUTable(listOfTaxa,distMatrix)
    deltaMatrix = produceDeltaMatrix(listOfTaxa,distMatrix,uTable)
    
    print()
    counter = 0 
    if N == 2 : 
        print ("Finally :")
        print(listOfTaxa)
    # find the smallest value in delta matrix 
    else :

        min_key = min(deltaMatrix,key = deltaMatrix.get)
        print ("-Pulling away  {} and {} ".format(min_key[0],min_key[1]))
        
        # calculate the distance between min_key and the parent node 
        parent_node = min_key[0] + min_key[1]
        dist1 = float(distMatrix[min_key]) + ((uTable[min_key[0]] - uTable[min_key[1]]) / (N-2) )
        dist1 = dist1/2

        print ("The distance between {} and  {}  is {} ".format(min_key[0],parent_node,dist1))
        dist2 = float(distMatrix[min_key]) - dist1

        print ("The distance between {} and  {}  is {} ".format(min_key[1],parent_node,dist2))

        # remove the taxas that were pulled out of the list of taxa's 
        listOfTaxa.remove(min_key[0])
        listOfTaxa.remove(min_key[1])
        # loop throgh the rest of taxa list and calculate the distance
        for values in listOfTaxa:
            # calculate the distance between the values and the parent node
            ik = (min_key[0], values)
            jk = (min_key[1],values)


            newdist = 1/2 * (float(distMatrix[ik] ) + float(distMatrix[jk]) - float(distMatrix[min_key]))

            result_key = (values, parent_node)
            distMatrix[result_key] = newdist
        


        # update  the distance matrix by removing all items od min_ky 
        for i in list(distMatrix):
            if (i[0] == min_key[0] or i[1] == min_key[0]) or (i[0] == min_key[1] or i[1] == min_key[1]) :
                del distMatrix[i]
        # add the new taxa 
        listOfTaxa.append(parent_node)

        #print(distMatrix)
        print()
        neighborJoiningRec(listOfTaxa, distMatrix)
        #print(distMatrix)


############
####MAIN####
############

if __name__ == '__main__':

    
    if len(sys.argv) != 2:
        print("USAGE:  python  NJ.py  distMatrix")
        sys.exit(0)

    taxaList, distDict = readDistMatrixFile(sys.argv[1])

    neighborJoiningRec(taxaList, distDict)

    # sm = produceUTable(taxaList,distDict)
    # tx = produceDeltaMatrix(taxaList,distDict,sm)
    # print (tx)
