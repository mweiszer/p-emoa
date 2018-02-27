#12.04.2016
#Michal Weiszer, University of Lincoln, UK
#Use as ParetoFilter(obj,n) where obj is a list of solutions, e.g. obj=[[3073.0, 875.0, 1984.0], [3073.0, 875.0, 1984.0], [3075.0, 874.0, 1984.0]] and n is the desired number of filtered solutions.
import scipy.io
from math import sqrt
import numpy
def distance(x,y):
    return sqrt(sum([(coord-y[i])**2 for i,coord in enumerate(x)]))

def midpoint(point1,point2):
    if len(point1)<>len(point2):
        raise ValueError('len of points is not the same')
    return [(point1[i]+point2[i])/2.0 for i in range(len(point1))]

def circles(point,solutions):



    max_diameter = 0
    min_diameter = float('inf')
    max_d_point = []
    min_d_point = []
    for sol1 in solutions:
        if sol1==point:
            continue

        diameter = distance(point,sol1)

        c = midpoint(sol1,point)
        for sol2 in solutions:
            if sol1==sol2:
                continue

            if distance(c,sol2)+0.0000001<(diameter/2):
                break
        else:
            if max_diameter<diameter:
                max_diameter = diameter
                max_d_point = sol1


            if min_diameter>diameter:
                min_diameter = diameter
                min_d_point = sol1
#                print "mind", sol1, diameter, c
    if max_diameter==0:
        print 'maxd=0', point

    return max_diameter, min_diameter

def diameters(sols):
    max_diameters = []
    min_diameters = []
    for p in sols:
        max_diameter = 0
        min_diameter = float('inf')
        max_d_point = []
        min_d_point = []

        for sol1 in sols:
            if sol1==p:
                continue

            diameter = distance(p,sol1)
            c = midpoint(sol1,p)
            for sol2 in sols:
                if sol1==sol2:
                    continue
                if distance(c,sol2)+0.0000001<(diameter/2):
                    break
            else:
                if max_diameter<diameter:
                    max_diameter = diameter
                    max_d_point = sol1

                if min_diameter>diameter:
                    min_diameter = diameter
                    min_d_point = sol1
        if max_diameter==0:
            print 'maxd=0', p
        max_diameters.append(max_diameter)
        min_diameters.append(min_diameter)

    avg_max_diameters = sum(max_diameters)/len(max_diameters)
    avg_min_diameters = sum(min_diameters)/len(min_diameters)
    return avg_max_diameters, avg_min_diameters
def xi_measure(sols):
    all_diameters = []
    
    for p in sols:
        max_diameter = 0
        min_diameter = float('inf')
        max_d_point = []
        min_d_point = []
        for sol1 in sols:
            if sol1==p:
                continue

            diameter = distance(p,sol1)
            c = midpoint(sol1,p)
            for sol2 in sols:
                if sol1==sol2:
                    continue
                if distance(c,sol2)+0.0000001<(diameter/2):
                    break
            else:
                if max_diameter<diameter:
                    max_diameter = diameter
                    max_d_point = sol1

                if min_diameter>diameter:
                    min_diameter = diameter
                    min_d_point = sol1
        if max_diameter==0:
            print 'maxd=0', p
        all_diameters.append(max_diameter)
        all_diameters.append(min_diameter)

 
    return numpy.std(all_diameters)/numpy.mean(all_diameters)
def dif_diameter(d):
    return (d[0]+d[1])/2
def even_tdea_archiver(population,n,threshold_n):

    

    new_archive = []
    pop_indexes = []
    max_diameters = []
    min_diameters = []
    new_archive.append(extreme_sol(population,'min',0))
    new_archive.append(extreme_sol(population,'min',1))
    pop_indexes.append(population.index(extreme_sol(population,'min',0)))
    pop_indexes.append(population.index(extreme_sol(population,'min',1)))
    for pop_index,ind in enumerate(population):
        
        if len(new_archive) == 0:
            
            new_archive.append(ind)
            pop_indexes.append(pop_index)
        else:
            should_remove = []
            should_add = True
            for a in new_archive:
                if ind == a:
                    should_add = False
                    break
                
            if should_add and len(new_archive) < 2:
                
                new_archive.append(ind)
                pop_indexes.append(pop_index)
            elif should_add and len(new_archive)<n:
                #tdea addition
                max_obj = [float('-inf') for _ in ind]
                min_obj = [float('inf') for _ in ind]
                a_star = []
                dist_min = float('inf')
                for i,a in enumerate(new_archive):
                    rect_dist = 0
                    for obj_i,obj in enumerate(a): 
                        rect_dist += abs(ind[obj_i]-obj)
                        if obj > max_obj[obj_i]:
                            max_obj[obj_i] = obj
                        if obj < min_obj[obj_i]:
                            min_obj[obj_i] = obj
                    if rect_dist<dist_min:
                        a_star = a
                        dist_min = rect_dist
                    
                delta = max([abs(ind[i]-obj) for i,obj in enumerate(a_star)])
                threshold = 1.0/(threshold_n)**(1.0/3)
                if delta>=threshold:
                    new_archive.append(ind)
                    pop_indexes.append(pop_index)
    if len(new_archive)<n:
        return [],new_archive
    if len(max_diameters)==0:
        for p in new_archive:
            max_diameter = 0
            min_diameter = float('inf')
            max_d_point = []
            min_d_point = []
            for sol1 in new_archive:
                if sol1==p:
                    continue

                diameter = distance(p,sol1)
                c = midpoint(sol1,p)
                for sol2 in new_archive:
                    if sol1==sol2:
                        continue
                    if distance(c,sol2)+0.0000001<(diameter/2):
                        break
                else:
                    if max_diameter<diameter:
                        max_diameter = diameter
                        max_d_point = sol1

                    if min_diameter>diameter:
                        min_diameter = diameter
                        min_d_point = sol1
            if max_diameter==0:
                print 'maxd=0', p
            max_diameters.append(max_diameter)
            min_diameters.append(min_diameter)


    avg_max_diameters = sum(max_diameters)/len(max_diameters)
    avg_min_diameters = sum(min_diameters)/len(min_diameters)
    avg_d = sum(max_diameters+min_diameters)/(len(max_diameters)*2)
    for pop_index,ind in enumerate(population):
        
        if len(new_archive) == 0:
            
            new_archive.append(ind)
            pop_indexes.append(i)
        else:
            should_remove = []
            should_add = True
            for a in new_archive:
                if ind == a:
                    should_add = False
                    break
            if should_add:
                #eveness measure

                nearest_dist = float('inf')
                nearest = []
                for sol1 in new_archive:
                    d=distance(ind,sol1)
                    if d<nearest_dist:
                        nearest_dist=d
                        nearest = sol1
                
                if abs(min(circles(ind,[s for s in new_archive if s <> nearest]))-avg_d) < abs(min(circles(nearest,[s for s in new_archive]))-avg_d) and abs(max(circles(ind,[s for s in new_archive if s <> nearest]))-avg_d) < abs(max(circles(nearest,[s for s in new_archive]))-avg_d):
                    new_archive.remove(nearest)
                    new_archive.append(ind)
                    
    return pop_indexes, new_archive

def normalisation(pop):
    """normalises the population for nay number of objectives"""
    pop_n = []
    maxs = [max([p[i] for p in pop]) for i in range(0,len(pop[0]))]
    mins = [min([p[i] for p in pop]) for i in range(0,len(pop[0]))]
    
    for sol in pop:
        sol_n = []
        for obj in range(0,len(pop[0])):
            sol_n.append((sol[obj]-mins[obj])/(maxs[obj]-mins[obj]))
        pop_n.append(sol_n)
    return pop_n
def denormalisation(pop_n,pop):
    """denormalises the population for nay number of objectives"""
    pop_dn = []
    maxs = [max([p[i] for p in pop]) for i in range(0,len(pop[0]))]
    mins = [min([p[i] for p in pop]) for i in range(0,len(pop[0]))]
    
    for sol in pop_n:
        sol_dn = []
        for obj in range(0,len(pop_n[0])):
            sol_dn.append(sol[obj]*(maxs[obj]-mins[obj])+mins[obj])
        pop_dn.append(sol_dn)
    return pop_dn
def normalisation_inc(pop, m=[], inc_n=[], tau=1):
    """normalises the population for any number of objectives with a fixed increment"""
    pop_n = []
    maxs = [max([p[i] for p in pop]) for i in range(0,len(pop[0]))]
    mins = [min([p[i] for p in pop]) for i in range(0,len(pop[0]))]
    
    for sol in pop:
        sol_n = []
        for obj in range(0,len(pop[0])):
            if obj not in m:
                adj = (inc_n[obj]/tau)-(maxs[obj]-mins[obj])
            elif obj in m:
                adj=0
            sol_n.append((sol[obj]-mins[obj])/(maxs[obj]-mins[obj]+adj))
        pop_n.append(sol_n)
    return pop_n  
def denormalisation_inc(pop_n, pop, m=[], inc_n=[],tau=1):
    """denormalises the population for any number of objectives with a fixed increment"""
    pop_dn = []
    maxs = [max([p[i] for p in pop]) for i in range(0,len(pop[0]))]
    mins = [min([p[i] for p in pop]) for i in range(0,len(pop[0]))]
    
    for sol in pop_n:
        sol_dn = []
        for obj in range(0,len(pop_n[0])):
            if obj not in m:
                adj = (inc_n[obj]/tau)-(maxs[obj]-mins[obj])
            elif obj in m:
                adj=0
            sol_dn.append((sol[obj])*(maxs[obj]-mins[obj]+adj)+mins[obj])
        pop_dn.append(sol_dn)
    return pop_dn  
def extreme_sol(sols,type,obj):
    sols0 = [s[obj] for s in sols]
    if type=='min':
        return sols[sols0.index(min(sols0))]
    elif type=='max':
        return sols[sols0.index(max(sols0))]

def ParetoFilter(obj,n):
    
    
    obj_n=normalisation(obj)
    
    
        
    i=0
    j = 0
    pop_indexes,sols=even_tdea_archiver(obj_n,n,n*(2**3))
    try:
        while len(sols)<n:
            i=i+1
            pop_indexes,sols=even_tdea_archiver(obj_n,n,n*(2**i))
            if len(sols)>=n:
        
                i = i-1
                pop_indexes,sols=even_tdea_archiver(obj_n,n,n*(2**i))
                while len(sols)<n:
                    i=i+0.001
                    pop_indexes,sols=even_tdea_archiver(obj_n,n,n*(2**i))
                i = i -0.001
                pop_indexes,sols=even_tdea_archiver(obj_n,n,n*(2**i))
                while len(sols)<n:
                    j=j+1
                    pop_indexes,sols=even_tdea_archiver(obj_n,n,n*(2**i)+j)
                break
            if i>100000:
                raise ValueError('Iterations exceeded')
        sols_dn = denormalisation(sols,obj)
    except:
        
        sols_dn = denormalisation(sols,obj)
    return sols_dn,pop_indexes