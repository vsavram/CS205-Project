import numpy as np
from mpi4py import MPI


def K_Means_parallel(data, comm, rank, mpi_size, k=2, tol=0.001, max_iter=300):


    centroids = {}

    for i in range(k):
        centroids[i] = data[i]

    for i in range(max_iter):
        classifications = {}

        for j in range(k):
            classifications[j] = []
 
        # scattering data (got from here: https://www.kth.se/blogs/pdc/2019/11/parallel-programming-in-python-mpi4py-part-2/)
        sendbuf = None
        if rank == 0:
            sendbuf = data
            ave, res = divmod(sendbuf.size, mpi_size)

            # count: the size of each sub-task
            count = [ave + 1 if p < res else ave for p in range(mpi_size)]
            count = np.array(count)

            # displacement: the starting index of each sub-task
            displ = [sum(count[:p]) for p in range(mpi_size)]
            displ = np.array(displ)
      
        else:
            sendbuf = None
            # initialize count on worker processes
            count = np.zeros(mpi_size, dtype=np.int)
            displ = None 

        # broadcast count
        comm.Bcast(count, root=0)

        # initialize recvbuf on all processes
        recvbuf = np.zeros(count[rank])

        comm.Scatterv([sendbuf, count, displ, MPI.DOUBLE], recvbuf, root=0)
        
        # broacast centroids
        centroids = comm.bcast(centroids, root = 0)
 
        #reshape data to array
        recvbuf = np.array(recvbuf).reshape(int(len(recvbuf)/data.shape[1]),data.shape[1])

        # loop each process does
        for featureset in recvbuf:
            distances = []
            for centroid in centroids:
                    distances.append(np.linalg.norm(featureset-centroids[centroid]))
                    #for ii,val in enumerate(featureset):
                    #    norm += (featureset[ii] - centroids[centroid][ii])**2
                    #distances.append(norm)
            classification = distances.index(min(distances))
            classifications[classification].append(featureset)
        
        #if i == max_iter-1:
        #    print('rank: ',rank,len(classifications[0]),len(classifications[1]))
        print(classifications.keys())
        classifications = comm.gather(classifications, root = 0)

        if rank == 0:
            
            #recreate classifications dict from gathered processes
            new_cls = {}
            for key in classifications[0].keys():
                print(list(d[key] for d in classifications))
                new_cls[key] = np.concatenate(list(d[key] for d in classifications),axis=0)
           
            classifications = new_cls 
            prev_centroids = dict(centroids)
    
            for classification in classifications:
                centroids[classification] = np.average(classifications[classification],axis=0)

            optimized = True

            for c in centroids:
                original_centroid = prev_centroids[c]
                current_centroid = centroids[c]
                if np.sum((current_centroid-original_centroid)/original_centroid*100.0) > tol:
                    #print(np.sum((current_centroid-original_centroid)/original_centroid*100.0))
                    optimized = False


    if rank == 0:
        if optimized:
            print("clustering is optimized")

        else:
            print("cluster is not optimized (possibly increase iterations?)")

    return centroids,classifications


if __name__ == "__main__":

    import pandas as pd
    import time
    import matplotlib.pyplot as plt

    # start MPI processes
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    mpi_size = comm.Get_size()
    
    exp_data = pd.read_csv("counts_small_test.csv",index_col=0)
    exp_data.drop(['Unnamed: 0.1'],axis=1,inplace=True)

    start_time = time.time()
    centroids, classifications = K_Means_parallel(exp_data[exp_data.columns[0:30]].values,
                                                        comm=comm,rank=rank,mpi_size=mpi_size,
                                                        k=2,max_iter=10)
    end_time = time.time()

    if rank == 0:
        print("\nExecution time: {}\n".format(end_time-start_time))
        print('Cluster sizes')
        print(['cluster {}: {}'.format(key,len(classifications[key])) for key in classifications.keys()])

        show_plots = False # only works when data dimensionality is 2
        if show_plots:
            for centroid in centroids:
                plt.scatter(centroids[centroid][0], centroids[centroid][1],
                            marker="o", color="k")

            for classification in classifications:
                for featureset in classifications[classification]:
                    plt.scatter(featureset[0], featureset[1], marker="x", color='b')

            plt.show()

