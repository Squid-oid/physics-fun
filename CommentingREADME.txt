Comments describing the nature of a function should universally take the following form,

def myfunc(param1:float = something, param2 = somethingelse, param3:custclass.custobj) -> sometype: 
    """
    Initial description of what the function does, can include surface level description of arguments or return

    Parameters
    ----------
    param1 : float
        A description of what param1 represents, and sometimes a hint of how to make one and a mention of that if none is given it will assume value something 
        which has some implication.
    param2 : type hint (ex. if param 2 has no type hint because it can be one of two types you could write type1/type2, or if it is completely free write FREE or ANY etc.)
        Same as for param1, perhaps include why it has no type hint
    param3 : custclass.custobj
        Mandatory, then write the description as in param1 and param2.
    
    Returns
    -------
    returnName : sometype
        Explain what is returned
    """

    ...body of myfunc...


An example from W_object.Ball during development

def bar_overlap(self:'Ball', bar:Barrier):
        '''
        Finds and returns the greatest overlapping distance between a ball and a barrier then returns it in the standard basis vectors as a vector

        Parameters
        ----------
        bar : W_object.Barrier
            The barrier to check overlap with
        Returns
        -------
        overlap : numpy.mat 
            The maximum overlap in vector form, directed such that a negative overlap implies that the two are not touching yet
        '''
        delta = CVector(self.coord.delta(bar.coord))
        proj_delta = CVector([delta.transform(bar.mat)[0,0], 0])
        proj_gap = proj_delta - CVector([self.radius, 0])
        overlap = np.matmul(proj_gap.vec, bar.mat)
        return overlap

An example for if you return a list

def flip_and_list(obj1, obj2) -> List:
    """
    Returns a list of the two input objects passed to the function

    Parameters
    ----------
    obj1 : ANY
        Any object, this will end up returned at the end of the list
    obj2 : ANY
        Any object, this will end up returned at the head of the list

    Returns
    -------
    returnList[a,b] : List
        A list of the input objects flipped
    a : ANY
        The second argument passed to the function
    b : ANY
        The first argument passed to the function
    """ 
    return_list = [obj2, obj1]
    return return_list

Of course the comments can be adapted to a situation for instance this is a function I wrote in a FEM course that has been slightly adapted to work mostly abide by this format
while maintaing legibility

def mesh_and_boundaries(geom:cfg.Geometry, id_qnh:int = 100, id_conv:int = 200, show_fig:bool = True, size_factor:float = 2.5):
    """
    Returns the nessecary matrices to begin work in the presolver. 
    These matrices are suitably named [coords, edof, dofs, bdofs, elementmarkers]
    Also returns the interesting (qn_h and convection) boundary elements in two arrays
    that are suitably named [ele_qnh, ele_conv]
    Essentialy a fancy wrapper for the Calfem mesh.create() function

    Parameters
    ----------
        geom : cfg.Geometry
            The geometry to mesh, reasonably the one created using base_shape()
        id_qnh : int
            The id used when creating the geom to track prescribed flow boundaries
        id_conv : int
            The id used when creating the geom to track convective flow boundaries
        show_fig : bool
            Wether or not to present the geometry before continuing
        size_factor : float
            The scaling factor that will be passed to the mesh command after being scaled
            to our base length

    Returns
    -------
        A list containing in order:
            mesh:
                The mesh created by cfm.GmshMesh
                Used for plotting/drawing
            coords:          
                Node coordinates
                [[n0_x, n0_y, n0_z],
                [   ...           ],
                [nn_x, nn_y, nn_z]]
            edof:            
                Element topology
                [[el0_dof1, ..., el0_dofn],
                [          ...          ],
                [eln_dof1, ..., eln_dofn]]
            dofs:            
                Node dofs
                [[n0_dof1, ..., n0_dofn],
                [         ...         ],
                [nn_dof1, ..., nn_dofn]]
            bdofs:
                Boundary dofs. Dictionary containing lists of dofs for
                each boundary marker. Dictionary key = marker id.
            elementmarkers: 
                List of integer markers. Row i contains the marker of
                element i. Markers are similar to boundary markers and
                can be used to identify in which region an element lies.
            ele_qnh:
                List of "elements" (two node pairs) that lie on the qn_h
                boundary, useful later for calculating the heat flow into
                these elements for the f matrix
            ele_conv
                List of "elements" (two node pairs) that lie on the conv
                boundary, useful later for calculating the convective 
                additions to the f and K matrices
    """
    L, t, a, b, c, d, h = define_shape_constants()

    mesh = cfm.GmshMesh(geom)
    mesh.elType = 2 # Type of element
    mesh.dofsPerNode = 1 # Degrees of freedom per node
    mesh.elSizeFactor = L*size_factor # Element size Factor

    mesh.returnBoundaryElements = True
    coords, edof, dofs, bdofs, elementmarkers, boundaryelements = mesh.create()

    ele_qnh = boundaryelements[id_qnh]      #qnh "elements" (node pairs)
    ele_conv = boundaryelements[id_conv]      #convective "elements" (node pairs)

    if show_fig == True :
        cfv.figure()
        cfv.drawMesh(coords=coords, edof=edof, dofs_per_node=mesh.dofsPerNode, el_type=mesh.elType, filled=True, title="Meshed Shape, size_factor = 5")
        cfv.show()

    return mesh, coords, edof, dofs, bdofs, elementmarkers, boundaryelements, ele_qnh, ele_conv

