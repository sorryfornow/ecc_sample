from collections import *



def getMinCost(val, t_nodes, t_from, t_to):
    # Write your code here
    val = [i%2 for i in val]
    if sum(val)==0:
        return 0
    # if val.count(1)%2:
    #     return -1

    # find the distance from leaf to root
    parent = dict()
    son = defaultdict(int)
    n = len(t_from)
    for i in range(n):
        parent[t_to[i]] = t_from[i]
        son[t_from[i]]+=1
    cost = 0
    # parents
    t_from = set(t_from)
    # sons
    t_to = set(t_to)
    leaves = t_to.difference(t_from)
    while len(leaves)>0:
        # find the min edge
        new_leaves = []
        for i in leaves:
            # find parent
            if i in parent:
                j = parent[i]
                # has to cost
                if val[i-1] != 0:
                    # delete the edge
                    val[j-1] = 1^val[j-1]
                    cost += 1
                son[j]-=1
                if son[j]==0:
                    new_leaves.append(j)
        leaves = new_leaves
    return cost

print(getMinCost([3,1,2],3,[1,1],[2,3]))


def studyTikcode(accounts, events):
    n = len(accounts)
    for eT,eN,eNew in events:
        if eT == 1:
            accounts = [eN if _ < eN else _ for _ in accounts]
        else:
            accounts[eN] = eNew
    return accounts