def deleteEdge(src, dest):
    sigma_old = [], distance_old = [], trackLost = [], pairsDone = []
    sinks = deleteUpdate(dest, src, src)
    sources = deleteUpdate(src, dest, dest)
    for s in sinks:
        deleteUpdate(src, dest, s)
    for s in sources:
        deleteUpdate(dest, src, s)
    increaseBetweenness()
    return


def deleteUpdate(src, dest, z):
    workSet = {src -> dest}
    visitedVertices = {src}
    affectedVertices = {}
    while workSet != {}:
        {x -> y} = Pop(workSet)
        alt = cost(x, y) + distance(y, z)
        if alt > distance(x, z):
            if<x, z> not in sigma_old:
                distance_old(x, z) = distance(x, z)
                sigma_old(x, z) = sigma(x, z)
                reduceBetweenness(x, z)
                # find new path and update sigma, P and D..
                for all the neighbors v of x:
                    choose the v with smallest distance(v, z) + cost(x, v)
                sigma(x, z) = sigma(v, z)
                P[x][z] = P[v][z]
                alt = cost(x, v) + distance(v, z)
                # sigma(x, z) = 0
                # clear P[x][z]
            if [x, z] in pairsDone:
                remove [x, z] from pairsDone
            distance(x, z) = alt
        if alt == distance(x, z) and distance(x, z) != infinity:
            if [x, z] not in pairsDone:
                if <x, z> not in sigma_old:
                    reduceBetweenness(x, z)
                # if sigma(x, z) != 0:
                    sigma_old(x, z) = sigma(x, z)
                    sigma(x, z) = sigma(x, z) - sigma(x, src) * 1 * sigma(dest, z)
                    if sigma(x, z) == 0:
                    #find new path, update sigma, P and D


            insert [x, z] to pairsDone
            insert x to affectedVertices
            for u in pred(x) ....:
                #the same as before..




def reduceBetweenness(a, z):
    # the same as before..
    return

def increaseBetweenness(a, z):
    # the same as before..
    return
