class UnionFind:
    """
    Define a union-find data structure. Also known as disjoint-set.
    The union-find structure groups n elements into a collection of disjoint sets.
    The basic operations consist of
       FIND  the unique set that contains the given element
       UNION of two sets

    The data structure mantains a collection S={S_i} of disjoint sets where
    each of them has an element that works as the representative.

    This particular implementation can be thought as a disjoint-set forest,
    where each element can be thought as a vertex and each disjoint set
    as a tree.
    """

    def __init__(self,mergeF=None):
        """
        Constructor of the union-find data structure

        Attributes:
            `V` : A dictionary of vertices
            `Vl`: A list containing the dictionary keys
            `parent`: A list of the parent node of the vertex
            `size`: List of numbers of vertices hanging from each vertex
            `mergeF`: (Optional) merge function of the data when merging two trees
            `data`: (Optional) Dictionary with additional information carried by each vertex
        """
        self.V = dict()
        self.Vl = []
        self.parent = []
        self.size = []
        self.mergeF = mergeF
        self.data = dict()

    def add(self,v,data=None):
        """
        Add a new element `v` to the union-find

        Parameters:
            `v`: Vertex to be inserted
            `data`: (Optional) information associated to the vertex
        """
        if v not in self.V:
            i = len(self.V)
            self.V[v] = i
            self.Vl.append(v)
            self.parent.append(i)
            self.size.append(1)
            self.data[v] = data

    def find_parent(self,v):
        """
        Find the root of the tree where v is located.
        Update the parent information for each vertex

        Parameters:
            `v`: Vertex
        """
        i = self.V[v]
        p = i
        while self.parent[p]!=p:
            p = self.parent[p]
        while i!=p:
            i,j = self.parent[i],i
            self.parent[j] = p
        return self.Vl[p]

    def find_size(self,v):
        """
        Returns the number of vertices hanging below v

        Parameters:
            `v`: vertex
        """
        return self.size[self.V[self.find_parent(v)]]

    # returns (new root,merged root)
    def merge(self,u,v):
        """
        The tree containing vertex u is merged INTO the tree containing vertex v.
        """

        su = self.find_size(u)
        sv = self.find_size(v)
        pu = self.parent[self.V[u]]
        pv = self.parent[self.V[v]]
        if pu==pv:
            return (self.Vl[pu],None)
        d = self.mergeF(self.data[self.Vl[pu]],self.data[self.Vl[pv]])
        if sv<=su:
            pu,pv,su,sv = pv,pu,sv,su
        self.parent[pv] = pu
        self.size[pu] = su+sv
        self.data[self.Vl[pu]] = d
        return (self.Vl[pu],self.Vl[pv])

    def getData(self,v):
        return self.data[self.find_parent(v)]

def mergeF(a,b):
    m,ei,ej = max((a['max'],a['elderi'],a['elderj']),(b['max'],b['elderi'],b['elderj']))
    return {'max':m,'elderi':ei, 'elderj':ej}

def persistence(f):
    """
    Computes the 0D persistence of a given filtration. Filtration values are
    passed in an array-like structure `f`. These values are later sorted from
    high to low and we compute the persistence of connected components.

    Whenever two connected components merge, we record the filter value difference
    and store the merge point in a list `pairs`.

    Parameters:
        `f`: array-like with filtration values.

    Returns:
        `pairs`: list of tuples (persistence, death, birth). A connected-component
        with infinite persistence is appended at the end.
    """
    fi = []
    for i in range(f.shape[0]):
        for j in range(f.shape[1]):
            for k in range(f.shape[2]):
                fi.append((f[i,j,k],i,j,k))
    
    fi = sorted(fi, reverse=True)
    
    uf = UnionFind(mergeF)
    pairs = []
    for v, i, j, k in fi:
        uf.add((i,j,k),{'max':v,'elderi':i, 'elderj':j, 'elderk':k})
        
        top = (i-1, j, k)
        bot = (i+1, j, k)
        lef = (i, j-1, k)
        rig = (i, j+1, k)
        fro = (i, j, k-1)
        bac = (i, j, k+1)
        
        topb = top in uf.V
        botb = bot in uf.V
        lefb = lef in uf.V
        rigb = rig in uf.V
        frob = fro in uf.V
        bacb = bac in uf.V
        

        if topb:
            topp = uf.find_parent(top)
            topd = uf.getData(top)
            
        if lefb:
            lefp = uf.find_parent(lef)
            lefd = uf.getData(lef)
            
        if botb:
            botp = uf.find_parent(bot)
            botd = uf.getData(bot)
            
        if rigb:
            rigp = uf.find_parent(rig)
            rigd = uf.getData(rig)
            
        if frob:
            frop = uf.find_parent(fro)
            frod = uf.getData(fro)
        
        if bacb:
            bacp = uf.find_parent(bac)
            bacd = uf.getData(bac)
            
            
        if topb:
            uf.merge(top, (i,j,k))
        if lefb:
            uf.merge(lef, (i,j,k))
        if frob:
            uf.merge(fro, (i,j,k))
        if botb:
            uf.merge((i,j,k), bot)
        if rigb:
            uf.merge((i,j,k), rig)
        if bacb:
            uf.merge((i,j,k), bac)
        
        if ~frob and ~bacb:
            
            if lefb and rigb and ~topb and ~botb:
                if lefp != rigp:
                    d,ci,cj,ck = min((lefd['max'], lefd['elderi'], lefd['elderj'], lefd['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
            
            elif ~lefb and ~rigb and topb and botb:
                if topp != botp:
                    d,ci,cj,ck = min((topd['max'], topd['elderi'], topd['elderj'],topd['elderk']), (botd['max'],botd['elderi'],botd['elderj'],botd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
            
            elif ~lefb and rigb and topb and ~botb:
                if rigp != topp:
                    d,ci,cj,ck = min((rigd['max'], rigd['elderi'], rigd['elderj'],rigd['elderk']), (topd['max'],topd['elderi'],topd['elderj'],topd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
            
            elif lefb and ~rigb and topb and ~botb:
                if lefp != topp:
                    d,ci,cj,ck = min((lefd['max'], lefd['elderi'], lefd['elderj'],lefd['elderk']), (topd['max'],topd['elderi'],topd['elderj'],topd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
            
            elif lefb and ~rigb and ~topb and botb:
                if lefp != botp:
                    d,ci,cj,ck = min((lefd['max'], lefd['elderi'], lefd['elderj'],lefd['elderk']), (botd['max'],botd['elderi'],botd['elderj'],botd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
            
            elif ~lefb and rigb and ~topb and botb:
                if rigp != botp:
                    d,ci,cj,ck = min((rigd['max'], rigd['elderi'], rigd['elderj'],rigd['elderk']), (botd['max'],botd['elderi'],botd['elderj'],botd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
            
            elif lefb and rigb and topb and ~botb:
                if (lefp != rigp) and (rigp != topp) and (topp != lefp):
                    d,ci,cj,ck = min((lefd['max'], lefd['elderi'], lefd['elderj'],lefd['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
            
                    d,ci,cj,ck = min((rigd['max'], rigd['elderi'], rigd['elderj'],rigd['elderk']), (topd['max'],topd['elderi'],topd['elderj'],topd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                
                elif (lefp != rigp) and (rigp != topp) and (topp == lefp):
                    d,ci,cj,ck = min((lefd['max'], lefd['elderi'], lefd['elderj'],lefd['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                
                elif (lefp != rigp) and (rigp == topp) and (topp != lefp):
                    d,ci,cj,ck = min((lefd['max'], lefd['elderi'], lefd['elderj'],lefd['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                
                elif (lefp == rigp) and (rigp != topp) and (topp != lefp):
                    d,ci,cj,ck = min((topd['max'], topd['elderi'], topd['elderj'],topd['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
            
            elif lefb and ~rigb and topb and botb:
                if (lefp != topp) and (topp != botp) and (botp != lefp):
                    d,ci,cj,ck = min((topd['max'], topd['elderi'], topd['elderj'],topd['elderk']), (botd['max'],botd['elderi'],botd['elderj'],botd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                    d,ci,cj,ck = min((botd['max'], botd['elderi'], botd['elderj'],botd['elderk']), (lefd['max'],lefd['elderi'],lefd['elderj'],lefd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (lefp != topp) and (topp != botp) and (botp == lefp):
                    d,ci,cj,ck = min((topd['max'], topd['elderi'], topd['elderj'],topd['elderk']), (botd['max'],botd['elderi'],botd['elderj'],botd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (lefp != topp) and (topp == botp) and (botp != lefp):
                    d,ci,cj,ck = min((topd['max'], topd['elderi'], topd['elderj'],topd['elderk']), (lefd['max'],lefd['elderi'],lefd['elderj'],lefd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (lefp == topp) and (topp != botp) and (botp != lefp):
                    d,ci,cj,ck = min((topd['max'], topd['elderi'], topd['elderj'],topd['elderk']), (botd['max'],botd['elderi'],botd['elderj'],botd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
            elif lefb and rigb and ~topb and botb:
                if (lefp != rigp) and (rigp != botp) and (botp != lefp):
                    d,ci,cj,ck = min((rigd['max'], rigd['elderi'], rigd['elderj'],rigd['elderk']), (botd['max'],botd['elderi'],botd['elderj'],botd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                    d,ci,cj,ck = min((botd['max'], botd['elderi'], botd['elderj'],botd['elderk']), (lefd['max'],lefd['elderi'],lefd['elderj'],lefd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (lefp != rigp) and (rigp != botp) and (botp == lefp):
                    d,ci,cj,ck = min((rigd['max'], rigd['elderi'], rigd['elderj'],rigd['elderk']), (botd['max'],botd['elderi'],botd['elderj'],botd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (lefp != rigp) and (rigp == botp) and (botp != lefp):
                    d,ci,cj,ck = min((rigd['max'], rigd['elderi'], rigd['elderj'],rigd['elderk']), (lefd['max'],lefd['elderi'],lefd['elderj'],lefd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (lefp == rigp) and (rigp != botp) and (botp != lefp):
                    d,ci,cj,ck = min((rigd['max'], rigd['elderi'], rigd['elderj'],rigd['elderk']), (botd['max'],botd['elderi'],botd['elderj'],botd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
            elif ~lefb and rigb and topb and botb:
                if (rigp != topp) and (topp != botp) and (botp != rigp):
                    d,ci,cj,ck = min((topd['max'], topd['elderi'], topd['elderj'],topd['elderk']), (botd['max'],botd['elderi'],botd['elderj'],botd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                    d,ci,cj,ck = min((botd['max'], botd['elderi'], botd['elderj'],botd['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (rigp != topp) and (topp != botp) and (botp == rigp):
                    d,ci,cj,ck = min((topd['max'], topd['elderi'], topd['elderj'],topd['elderk']), (botd['max'],botd['elderi'],botd['elderj'],botd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (rigp != topp) and (topp == botp) and (botp != rigp):
                    d,ci,cj,ck = min((topd['max'], topd['elderi'], topd['elderj'],topd['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (rigp == topp) and (topp != botp) and (botp != rigp):
                    d,ci,cj,ck = min((topd['max'], topd['elderi'], topd['elderj'],topd['elderk']), (botd['max'],botd['elderi'],botd['elderj'],botd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
            elif lefb and rigb and topb and botb:
                if (lefp != rigp) and (rigp != topp) and (topp != botp) and (botp != lefp) and (lefp != topp) and (rigp != botp):
                    d,ci,cj,ck = min((lefd['max'], lefd['elderi'], lefd['elderj'],lefd['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                    d,ci,cj,ck = min((topd['max'], topd['elderi'], topd['elderj'],topd['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                    d,ci,cj,ck = min((botd['max'], botd['elderi'], botd['elderj'],botd['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (lefb != rigb) and (rigp != topp) and (topp != botp) and (botp != lefp) and (left != topp) and (rigp == botp):
                    d,ci,cj,ck = min((lefd['max'], lefd['elderi'], lefd['elderj'],lefd['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                    d,ci,cj,ck = min((topd['max'], topd['elderi'], topd['elderj'],topd['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                
                elif (lefb != rigb) and (rigp != topp) and (topp != botp) and (botp == lefp) and (left != topp) and (rigp != botp):
                    d,ci,cj,ck = min((lefd['max'], lefd['elderi'], lefd['elderj'],lefd['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                    d,ci,cj,ck = min((topd['max'], topd['elderi'], topd['elderj'],topd['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (lefb != rigb) and (rigp != topp) and (topp == botp) and (botp != lefp) and (left != topp) and (rigp != botp):
                    d,ci,cj,ck = min((lefd['max'], lefd['elderi'], lefd['elderj'],lefd['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                    d,ci,cj,ck = min((topd['max'], topd['elderi'], topd['elderj'],topd['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (lefb != rigb) and (rigp != topp) and (topp != botp) and (botp != lefp) and (left == topp) and (rigp != botp):
                    d,ci,cj,ck = min((lefd['max'], lefd['elderi'], lefd['elderj'],lefd['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                    d,ci,cj,ck = min((botd['max'], botd['elderi'], botd['elderj'],botd['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (lefb != rigb) and (rigp == topp) and (topp != botp) and (botp != lefp) and (left != topp) and (rigp != botp):
                    d,ci,cj,ck = min((lefd['max'], lefd['elderi'], lefd['elderj'],lefd['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                    d,ci,cj,ck = min((botd['max'], botd['elderi'], botd['elderj'],botd['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (lefb == rigb) and (rigp != topp) and (topp != botp) and (botp != lefp) and (left != topp) and (rigp != botp):
                    d,ci,cj,ck = min((lefd['max'], lefd['elderi'], lefd['elderj'],lefd['elderk']), (topd['max'],topd['elderi'],topd['elderj'],topd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                    d,ci,cj,ck = min((lefd['max'], lefd['elderi'], lefd['elderj'],lefd['elderk']), (botd['max'],botd['elderi'],botd['elderj'],botd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                
                elif (lefb != rigb) and (rigp == topp) and (topp != botp) and (botp == lefp) and (left != topp) and (rigp != botp):
                    d,ci,cj,ck = min((lefd['max'], lefd['elderi'], lefd['elderj'],lefd['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (lefb != rigb) and (rigp != topp) and (topp != botp) and (botp != lefp) and (left == topp) and (rigp == botp):
                    d,ci,cj,ck = min((lefd['max'], lefd['elderi'], lefd['elderj'],lefd['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (lefb != rigb) and (rigp != topp) and (topp == botp) and (botp == lefp) and (left == topp) and (rigp != botp):
                    d,ci,cj,ck = min((lefd['max'], lefd['elderi'], lefd['elderj'],lefd['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (lefb != rigb) and (rigp == topp) and (topp == botp) and (botp != lefp) and (left != topp) and (rigp == botp):
                    d,ci,cj,ck = min((lefd['max'], lefd['elderi'], lefd['elderj'],lefd['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                
                elif (lefb == rigb) and (rigp == topp) and (topp != botp) and (botp != lefp) and (left == topp) and (rigp != botp):
                    d,ci,cj,ck = min((topd['max'], topd['elderi'], topd['elderj'],topd['elderk']), (botd['max'],botd['elderi'],botd['elderj'],botd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (lefb == rigb) and (rigp != topp) and (topp != botp) and (botp == lefp) and (left != topp) and (rigp == botp):
                    d,ci,cj,ck = min((topd['max'], topd['elderi'], topd['elderj'],topd['elderk']), (botd['max'],botd['elderi'],botd['elderj'],botd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
        
        elif ~lefb and ~rigb:
            
            if frob and bacb and ~topb and ~botb:
                if frop != bacp:
                    d,ci,cj,ck = min((frod['max'], frod['elderi'], frod['elderj'], frod['elderk']), (bacd['max'],bacd['elderi'],bacd['elderj'],bacd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
            
            elif ~frob and ~bacb and topb and botb:
                if topp != botp:
                    d,ci,cj,ck = min((topd['max'], topd['elderi'], topd['elderj'],topd['elderk']), (botd['max'],botd['elderi'],botd['elderj'],botd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
            
            elif ~frob and bacb and topb and ~botb:
                if bacp != topp:
                    d,ci,cj,ck = min((bacd['max'], bacd['elderi'], bacd['elderj'],bacd['elderk']), (topd['max'],topd['elderi'],topd['elderj'],topd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
            
            elif frob and ~bacb and topb and ~botb:
                if frop != topp:
                    d,ci,cj,ck = min((frod['max'], frod['elderi'], frod['elderj'],frod['elderk']), (topd['max'],topd['elderi'],topd['elderj'],topd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
            
            elif frob and ~bacb and ~topb and botb:
                if frop != botp:
                    d,ci,cj,ck = min((frod['max'], frod['elderi'], frod['elderj'],frod['elderk']), (botd['max'],botd['elderi'],botd['elderj'],botd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
            
            elif ~frob and bacb and ~topb and botb:
                if bacp != botp:
                    d,ci,cj,ck = min((bacd['max'], bacd['elderi'], bacd['elderj'],bacd['elderk']), (botd['max'],botd['elderi'],botd['elderj'],botd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
            
            elif frob and bacb and topb and ~botb:
                if (frop != bacp) and (bacp != topp) and (topp != frop):
                    d,ci,cj,ck = min((frod['max'], frod['elderi'], frod['elderj'],frod['elderk']), (bacd['max'],bacd['elderi'],bacd['elderj'],bacd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
            
                    d,ci,cj,ck = min((bacd['max'], bacd['elderi'], bacd['elderj'],bacd['elderk']), (topd['max'],topd['elderi'],topd['elderj'],topd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                
                elif (frop != bacp) and (bacp != topp) and (topp == frop):
                    d,ci,cj,ck = min((frod['max'], frod['elderi'], frod['elderj'],frod['elderk']), (bacd['max'],bacd['elderi'],bacd['elderj'],bacd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                
                elif (frop != bacp) and (bacp == topp) and (topp != frop):
                    d,ci,cj,ck = min((frod['max'], frod['elderi'], frod['elderj'],frod['elderk']), (bacd['max'],bacd['elderi'],bacd['elderj'],bacd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                
                elif (frop == bacp) and (bacp != topp) and (topp != frop):
                    d,ci,cj,ck = min((topd['max'], topd['elderi'], topd['elderj'],topd['elderk']), (bacd['max'],bacd['elderi'],bacd['elderj'],bacd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
            
            elif frob and ~bacb and topb and botb:
                if (frop != topp) and (topp != botp) and (botp != frop):
                    d,ci,cj,ck = min((topd['max'], topd['elderi'], topd['elderj'],topd['elderk']), (botd['max'],botd['elderi'],botd['elderj'],botd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                    d,ci,cj,ck = min((botd['max'], botd['elderi'], botd['elderj'],botd['elderk']), (frod['max'],frod['elderi'],frod['elderj'],frod['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (frop != topp) and (topp != botp) and (botp == frop):
                    d,ci,cj,ck = min((topd['max'], topd['elderi'], topd['elderj'],topd['elderk']), (botd['max'],botd['elderi'],botd['elderj'],botd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (frop != topp) and (topp == botp) and (botp != frop):
                    d,ci,cj,ck = min((topd['max'], topd['elderi'], topd['elderj'],topd['elderk']), (frod['max'],frod['elderi'],frod['elderj'],frod['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (frop == topp) and (topp != botp) and (botp != frop):
                    d,ci,cj,ck = min((topd['max'], topd['elderi'], topd['elderj'],topd['elderk']), (botd['max'],botd['elderi'],botd['elderj'],botd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
            elif frob and bacb and ~topb and botb:
                if (frop != bacp) and (bacp != botp) and (botp != frop):
                    d,ci,cj,ck = min((bacd['max'], bacd['elderi'], bacd['elderj'],bacd['elderk']), (botd['max'],botd['elderi'],botd['elderj'],botd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                    d,ci,cj,ck = min((botd['max'], botd['elderi'], botd['elderj'],botd['elderk']), (frod['max'],frod['elderi'],frod['elderj'],frod['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (frop != bacp) and (bacp != botp) and (botp == frop):
                    d,ci,cj,ck = min((bacd['max'], bacd['elderi'], bacd['elderj'],bacd['elderk']), (botd['max'],botd['elderi'],botd['elderj'],botd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (frop != bacp) and (bacp == botp) and (botp != frop):
                    d,ci,cj,ck = min((bacd['max'], bacd['elderi'], bacd['elderj'],bacd['elderk']), (frod['max'],frod['elderi'],frod['elderj'],frod['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (frop == bacp) and (bacp != botp) and (botp != frop):
                    d,ci,cj,ck = min((bacd['max'], bacd['elderi'], bacd['elderj'],bacd['elderk']), (botd['max'],botd['elderi'],botd['elderj'],botd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
            elif ~frob and bacb and topb and botb:
                if (bacp != topp) and (topp != botp) and (botp != bacp):
                    d,ci,cj,ck = min((topd['max'], topd['elderi'], topd['elderj'],topd['elderk']), (botd['max'],botd['elderi'],botd['elderj'],botd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                    d,ci,cj,ck = min((botd['max'], botd['elderi'], botd['elderj'],botd['elderk']), (bacd['max'],bacd['elderi'],bacd['elderj'],bacd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (bacp != topp) and (topp != botp) and (botp == bacp):
                    d,ci,cj,ck = min((topd['max'], topd['elderi'], topd['elderj'],topd['elderk']), (botd['max'],botd['elderi'],botd['elderj'],botd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (bacp != topp) and (topp == botp) and (botp != bacp):
                    d,ci,cj,ck = min((topd['max'], topd['elderi'], topd['elderj'],topd['elderk']), (bacd['max'],bacd['elderi'],bacd['elderj'],bacd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (bacp == topp) and (topp != botp) and (botp != bacp):
                    d,ci,cj,ck = min((topd['max'], topd['elderi'], topd['elderj'],topd['elderk']), (botd['max'],botd['elderi'],botd['elderj'],botd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
            elif frob and bacb and topb and botb:
                if (frob != bacb) and (bacp != topp) and (topp != botp) and (botp != frop) and (frot != topp) and (bacp != botp):
                    d,ci,cj,ck = min((frod['max'], frod['elderi'], frod['elderj'],frod['elderk']), (bacd['max'],bacd['elderi'],bacd['elderj'],bacd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                    d,ci,cj,ck = min((topd['max'], topd['elderi'], topd['elderj'],topd['elderk']), (bacd['max'],bacd['elderi'],bacd['elderj'],bacd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                    d,ci,cj,ck = min((botd['max'], botd['elderi'], botd['elderj'],botd['elderk']), (bacd['max'],bacd['elderi'],bacd['elderj'],bacd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (frob != bacb) and (bacp != topp) and (topp != botp) and (botp != frop) and (frot != topp) and (bacp == botp):
                    d,ci,cj,ck = min((frod['max'], frod['elderi'], frod['elderj'],frod['elderk']), (bacd['max'],bacd['elderi'],bacd['elderj'],bacd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                    d,ci,cj,ck = min((topd['max'], topd['elderi'], topd['elderj'],topd['elderk']), (bacd['max'],bacd['elderi'],bacd['elderj'],bacd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                
                elif (frob != bacb) and (bacp != topp) and (topp != botp) and (botp == frop) and (frot != topp) and (bacp != botp):
                    d,ci,cj,ck = min((frod['max'], frod['elderi'], frod['elderj'],frod['elderk']), (bacd['max'],bacd['elderi'],bacd['elderj'],bacd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                    d,ci,cj,ck = min((topd['max'], topd['elderi'], topd['elderj'],topd['elderk']), (bacd['max'],bacd['elderi'],bacd['elderj'],bacd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (frob != bacb) and (bacp != topp) and (topp == botp) and (botp != frop) and (frot != topp) and (bacp != botp):
                    d,ci,cj,ck = min((frod['max'], frod['elderi'], frod['elderj'],frod['elderk']), (bacd['max'],bacd['elderi'],bacd['elderj'],bacd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                    d,ci,cj,ck = min((topd['max'], topd['elderi'], topd['elderj'],topd['elderk']), (bacd['max'],bacd['elderi'],bacd['elderj'],bacd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (frob != bacb) and (bacp != topp) and (topp != botp) and (botp != frop) and (frot == topp) and (bacp != botp):
                    d,ci,cj,ck = min((frod['max'], frod['elderi'], frod['elderj'],frod['elderk']), (bacd['max'],bacd['elderi'],bacd['elderj'],bacd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                    d,ci,cj,ck = min((botd['max'], botd['elderi'], botd['elderj'],botd['elderk']), (bacd['max'],bacd['elderi'],bacd['elderj'],bacd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (frob != bacb) and (bacp == topp) and (topp != botp) and (botp != frop) and (frot != topp) and (bacp != botp):
                    d,ci,cj,ck = min((frod['max'], frod['elderi'], frod['elderj'],frod['elderk']), (bacd['max'],bacd['elderi'],bacd['elderj'],bacd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                    d,ci,cj,ck = min((botd['max'], botd['elderi'], botd['elderj'],botd['elderk']), (bacd['max'],bacd['elderi'],bacd['elderj'],bacd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (frob == bacb) and (bacp != topp) and (topp != botp) and (botp != frop) and (frot != topp) and (bacp != botp):
                    d,ci,cj,ck = min((frod['max'], frod['elderi'], frod['elderj'],frod['elderk']), (topd['max'],topd['elderi'],topd['elderj'],topd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                    d,ci,cj,ck = min((frod['max'], frod['elderi'], frod['elderj'],frod['elderk']), (botd['max'],botd['elderi'],botd['elderj'],botd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                
                elif (frob != bacb) and (bacp == topp) and (topp != botp) and (botp == frop) and (frot != topp) and (bacp != botp):
                    d,ci,cj,ck = min((frod['max'], frod['elderi'], frod['elderj'],frod['elderk']), (bacd['max'],bacd['elderi'],bacd['elderj'],bacd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (frob != bacb) and (bacp != topp) and (topp != botp) and (botp != frop) and (frot == topp) and (bacp == botp):
                    d,ci,cj,ck = min((frod['max'], frod['elderi'], frod['elderj'],frod['elderk']), (bacd['max'],bacd['elderi'],bacd['elderj'],bacd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (frob != bacb) and (bacp != topp) and (topp == botp) and (botp == frop) and (frot == topp) and (bacp != botp):
                    d,ci,cj,ck = min((frod['max'], frod['elderi'], frod['elderj'],frod['elderk']), (bacd['max'],bacd['elderi'],bacd['elderj'],bacd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (frob != bacb) and (bacp == topp) and (topp == botp) and (botp != frop) and (frot != topp) and (bacp == botp):
                    d,ci,cj,ck = min((frod['max'], frod['elderi'], frod['elderj'],frod['elderk']), (bacd['max'],bacd['elderi'],bacd['elderj'],bacd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                
                elif (frob == bacb) and (bacp == topp) and (topp != botp) and (botp != frop) and (frot == topp) and (bacp != botp):
                    d,ci,cj,ck = min((topd['max'], topd['elderi'], topd['elderj'],topd['elderk']), (botd['max'],botd['elderi'],botd['elderj'],botd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (frob == bacb) and (bacp != topp) and (topp != botp) and (botp == frop) and (frot != topp) and (bacp == botp):
                    d,ci,cj,ck = min((topd['max'], topd['elderi'], topd['elderj'],topd['elderk']), (botd['max'],botd['elderi'],botd['elderj'],botd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
        
        elif ~topb and ~botb:
            
            if lefb and rigb and ~frob and ~bacb:
                if lefp != rigp:
                    d,ci,cj,ck = min((lefd['max'], lefd['elderi'], lefd['elderj'], lefd['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
            
            elif ~lefb and ~rigb and frob and bacb:
                if frop != bacp:
                    d,ci,cj,ck = min((frod['max'], frod['elderi'], frod['elderj'],frod['elderk']), (bacd['max'],bacd['elderi'],bacd['elderj'],bacd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
            
            elif ~lefb and rigb and frob and ~bacb:
                if rigp != frop:
                    d,ci,cj,ck = min((rigd['max'], rigd['elderi'], rigd['elderj'],rigd['elderk']), (frod['max'],frod['elderi'],frod['elderj'],frod['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
            
            elif lefb and ~rigb and frob and ~bacb:
                if lefp != frop:
                    d,ci,cj,ck = min((lefd['max'], lefd['elderi'], lefd['elderj'],lefd['elderk']), (frod['max'],frod['elderi'],frod['elderj'],frod['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
            
            elif lefb and ~rigb and ~frob and bacb:
                if lefp != bacp:
                    d,ci,cj,ck = min((lefd['max'], lefd['elderi'], lefd['elderj'],lefd['elderk']), (bacd['max'],bacd['elderi'],bacd['elderj'],bacd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
            
            elif ~lefb and rigb and ~frob and bacb:
                if rigp != bacp:
                    d,ci,cj,ck = min((rigd['max'], rigd['elderi'], rigd['elderj'],rigd['elderk']), (bacd['max'],bacd['elderi'],bacd['elderj'],bacd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
            
            elif lefb and rigb and frob and ~bacb:
                if (lefp != rigp) and (rigp != frop) and (frop != lefp):
                    d,ci,cj,ck = min((lefd['max'], lefd['elderi'], lefd['elderj'],lefd['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
            
                    d,ci,cj,ck = min((rigd['max'], rigd['elderi'], rigd['elderj'],rigd['elderk']), (frod['max'],frod['elderi'],frod['elderj'],frod['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                
                elif (lefp != rigp) and (rigp != frop) and (frop == lefp):
                    d,ci,cj,ck = min((lefd['max'], lefd['elderi'], lefd['elderj'],lefd['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                
                elif (lefp != rigp) and (rigp == frop) and (frop != lefp):
                    d,ci,cj,ck = min((lefd['max'], lefd['elderi'], lefd['elderj'],lefd['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                
                elif (lefp == rigp) and (rigp != frop) and (frop != lefp):
                    d,ci,cj,ck = min((frod['max'], frod['elderi'], frod['elderj'],frod['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
            
            elif lefb and ~rigb and frob and bacb:
                if (lefp != frop) and (frop != bacp) and (bacp != lefp):
                    d,ci,cj,ck = min((frod['max'], frod['elderi'], frod['elderj'],frod['elderk']), (bacd['max'],bacd['elderi'],bacd['elderj'],bacd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                    d,ci,cj,ck = min((bacd['max'], bacd['elderi'], bacd['elderj'],bacd['elderk']), (lefd['max'],lefd['elderi'],lefd['elderj'],lefd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (lefp != frop) and (frop != bacp) and (bacp == lefp):
                    d,ci,cj,ck = min((frod['max'], frod['elderi'], frod['elderj'],frod['elderk']), (bacd['max'],bacd['elderi'],bacd['elderj'],bacd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (lefp != frop) and (frop == bacp) and (bacp != lefp):
                    d,ci,cj,ck = min((frod['max'], frod['elderi'], frod['elderj'],frod['elderk']), (lefd['max'],lefd['elderi'],lefd['elderj'],lefd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (lefp == frop) and (frop != bacp) and (bacp != lefp):
                    d,ci,cj,ck = min((frod['max'], frod['elderi'], frod['elderj'],frod['elderk']), (bacd['max'],bacd['elderi'],bacd['elderj'],bacd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
            elif lefb and rigb and ~frob and bacb:
                if (lefp != rigp) and (rigp != bacp) and (bacp != lefp):
                    d,ci,cj,ck = min((rigd['max'], rigd['elderi'], rigd['elderj'],rigd['elderk']), (bacd['max'],bacd['elderi'],bacd['elderj'],bacd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                    d,ci,cj,ck = min((bacd['max'], bacd['elderi'], bacd['elderj'],bacd['elderk']), (lefd['max'],lefd['elderi'],lefd['elderj'],lefd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (lefp != rigp) and (rigp != bacp) and (bacp == lefp):
                    d,ci,cj,ck = min((rigd['max'], rigd['elderi'], rigd['elderj'],rigd['elderk']), (bacd['max'],bacd['elderi'],bacd['elderj'],bacd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (lefp != rigp) and (rigp == bacp) and (bacp != lefp):
                    d,ci,cj,ck = min((rigd['max'], rigd['elderi'], rigd['elderj'],rigd['elderk']), (lefd['max'],lefd['elderi'],lefd['elderj'],lefd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (lefp == rigp) and (rigp != bacp) and (bacp != lefp):
                    d,ci,cj,ck = min((rigd['max'], rigd['elderi'], rigd['elderj'],rigd['elderk']), (bacd['max'],bacd['elderi'],bacd['elderj'],bacd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
            elif ~lefb and rigb and frob and bacb:
                if (rigp != frop) and (frop != bacp) and (bacp != rigp):
                    d,ci,cj,ck = min((frod['max'], frod['elderi'], frod['elderj'],frod['elderk']), (bacd['max'],bacd['elderi'],bacd['elderj'],bacd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                    d,ci,cj,ck = min((bacd['max'], bacd['elderi'], bacd['elderj'],bacd['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (rigp != frop) and (frop != bacp) and (bacp == rigp):
                    d,ci,cj,ck = min((frod['max'], frod['elderi'], frod['elderj'],frod['elderk']), (bacd['max'],bacd['elderi'],bacd['elderj'],bacd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (rigp != frop) and (frop == bacp) and (bacp != rigp):
                    d,ci,cj,ck = min((frod['max'], frod['elderi'], frod['elderj'],frod['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (rigp == frop) and (frop != bacp) and (bacp != rigp):
                    d,ci,cj,ck = min((frod['max'], frod['elderi'], frod['elderj'],frod['elderk']), (bacd['max'],bacd['elderi'],bacd['elderj'],bacd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
            elif lefb and rigb and frob and bacb:
                if (lefb != rigb) and (rigp != frop) and (frop != bacp) and (bacp != lefp) and (left != frop) and (rigp != bacp):
                    d,ci,cj,ck = min((lefd['max'], lefd['elderi'], lefd['elderj'],lefd['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                    d,ci,cj,ck = min((frod['max'], frod['elderi'], frod['elderj'],frod['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                    d,ci,cj,ck = min((bacd['max'], bacd['elderi'], bacd['elderj'],bacd['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (lefb != rigb) and (rigp != frop) and (frop != bacp) and (bacp != lefp) and (left != frop) and (rigp == bacp):
                    d,ci,cj,ck = min((lefd['max'], lefd['elderi'], lefd['elderj'],lefd['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                    d,ci,cj,ck = min((frod['max'], frod['elderi'], frod['elderj'],frod['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                
                elif (lefb != rigb) and (rigp != frop) and (frop != bacp) and (bacp == lefp) and (left != frop) and (rigp != bacp):
                    d,ci,cj,ck = min((lefd['max'], lefd['elderi'], lefd['elderj'],lefd['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                    d,ci,cj,ck = min((frod['max'], frod['elderi'], frod['elderj'],frod['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (lefb != rigb) and (rigp != frop) and (frop == bacp) and (bacp != lefp) and (left != frop) and (rigp != bacp):
                    d,ci,cj,ck = min((lefd['max'], lefd['elderi'], lefd['elderj'],lefd['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                    d,ci,cj,ck = min((frod['max'], frod['elderi'], frod['elderj'],frod['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (lefb != rigb) and (rigp != frop) and (frop != bacp) and (bacp != lefp) and (left == frop) and (rigp != bacp):
                    d,ci,cj,ck = min((lefd['max'], lefd['elderi'], lefd['elderj'],lefd['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                    d,ci,cj,ck = min((bacd['max'], bacd['elderi'], bacd['elderj'],bacd['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (lefb != rigb) and (rigp == frop) and (frop != bacp) and (bacp != lefp) and (left != frop) and (rigp != bacp):
                    d,ci,cj,ck = min((lefd['max'], lefd['elderi'], lefd['elderj'],lefd['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                    d,ci,cj,ck = min((bacd['max'], bacd['elderi'], bacd['elderj'],bacd['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (lefb == rigb) and (rigp != frop) and (frop != bacp) and (bacp != lefp) and (left != frop) and (rigp != bacp):
                    d,ci,cj,ck = min((lefd['max'], lefd['elderi'], lefd['elderj'],lefd['elderk']), (frod['max'],frod['elderi'],frod['elderj'],frod['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                    d,ci,cj,ck = min((lefd['max'], lefd['elderi'], lefd['elderj'],lefd['elderk']), (bacd['max'],bacd['elderi'],bacd['elderj'],bacd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                
                elif (lefb != rigb) and (rigp == frop) and (frop != bacp) and (bacp == lefp) and (left != frop) and (rigp != bacp):
                    d,ci,cj,ck = min((lefd['max'], lefd['elderi'], lefd['elderj'],lefd['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (lefb != rigb) and (rigp != frop) and (frop != bacp) and (bacp != lefp) and (left == frop) and (rigp == bacp):
                    d,ci,cj,ck = min((lefd['max'], lefd['elderi'], lefd['elderj'],lefd['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (lefb != rigb) and (rigp != frop) and (frop == bacp) and (bacp == lefp) and (left == frop) and (rigp != bacp):
                    d,ci,cj,ck = min((lefd['max'], lefd['elderi'], lefd['elderj'],lefd['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (lefb != rigb) and (rigp == frop) and (frop == bacp) and (bacp != lefp) and (left != frop) and (rigp == bacp):
                    d,ci,cj,ck = min((lefd['max'], lefd['elderi'], lefd['elderj'],lefd['elderk']), (rigd['max'],rigd['elderi'],rigd['elderj'],rigd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                
                elif (lefb == rigb) and (rigp == frop) and (frop != bacp) and (bacp != lefp) and (left == frop) and (rigp != bacp):
                    d,ci,cj,ck = min((frod['max'], frod['elderi'], frod['elderj'],frod['elderk']), (bacd['max'],bacd['elderi'],bacd['elderj'],bacd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
                    
                elif (lefb == rigb) and (rigp != frop) and (frop != bacp) and (bacp == lefp) and (left != frop) and (rigp == bacp):
                    d,ci,cj,ck = min((frod['max'], frod['elderi'], frod['elderj'],frod['elderk']), (bacd['max'],bacd['elderi'],bacd['elderj'],bacd['elderk']))
                    pairs.append((d-v, (i,j,k), (ci,cj,ck)))
        
        elif
                
    pairs.append((float('inf'),None,(fi[0][1], fi[0][2])))
    
    return pairs

def rel_persistence(f,threshold=1e4):
    """
    Same idea as `persistence`, except that instead of recording persistance
    (death - birth), we record persistence relative to death, that is
    (death - birth)/death.

    We can also limit ourselves to only record those critical points where
    persistence is larger than some set threshold in order to avoid noise such
    as birth=0, death=3, which results in a relative persistence of 1.

    Parameters:
        `f`: array-like with filtration values.
        `threshold`: scalar. Consider only critical points whose persistence is larger
        than a fixed threshold.

    Returns:
        `pairs`: list of tuples (persistence/death, death, birth).
        A connected-component with infinite persistence is appended at the end.
    """

    fi = sorted(list(zip(f,range(len(f)))),reverse=True)
    uf = UnionFind(mergeF)
    pairs = []
    for v,i in fi:
        uf.add(i,{'max':v,'elder':i})
        if i-1 in uf.V and i+1 in uf.V:
            a = uf.getData(i-1)
            b = uf.getData(i+1)
            d,j = min((a['max'],a['elder']),(b['max'],b['elder']))
            if d-v > threshold:
                pairs.append(((d-v)/d,i,j))
        if i-1 in uf.V:
            uf.merge(i-1,i)
        if i+1 in uf.V:
            uf.merge(i,i+1)
    pairs.append((float('inf'),None,fi[0][1]))
    return pairs
