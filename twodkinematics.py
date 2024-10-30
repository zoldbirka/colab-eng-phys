

# Approximate derivative calculation

def deriv(xx_tab:list, ff_tab:list) -> list:  # general math formalism: the function is f(x)
    '''
    Usage:
    input: xx_tab, ff_tab : two list of floats with the same length
            xx_tab: list of values, usually with a constant delta_x step size, 
                    example: xx_tab=[0.0, 0.1, 0.2, 0.3, ...]
            ff_tab: values of a function f(x) at the values of xx_tab
                    example: ff_tab=[f(0.0), f(0.1), f(0.2), ...]
    output: dff_xx: approximate derivative of f(x) at the points of xx_tab:
                    dff_xx=[df/dx(x=0.0), df/dx(x=0.1), df/dx(x=0.2), ...]

    Example usage: calculating velocity from position data.
                    vx_tab=deriv(t_tab, x_tab)
    '''
    N = len(xx_tab)    # size of array
    dff_dxx = [0.0]*N  # array of derivatives
    
    dff_dxx[0] = (ff_tab[1]-ff_tab[0])/(xx_tab[1]-xx_tab[0])  # first element
    dff_dxx[N-1] = (ff_tab[N-1]-ff_tab[N-2])/(xx_tab[N-1]-xx_tab[N-2]) # last element
        
    for i in range(1, N-1):  # middle elements
        dff_dxx[i] = (ff_tab[i+1]-ff_tab[i-1])/(xx_tab[i+1]-xx_tab[i-1])
        
    # corrections at the edges
    dff_dxx[0] = 2*dff_dxx[0] - dff_dxx[1]  
    dff_dxx[N-1] = 2*dff_dxx[N-1] - dff_dxx[N-2]
    return(dff_dxx)


# Approximate integral calculation

def integral(xx_tab:list, ff_tab:list, F0 = 0.0)->list:  # F0: start value
    '''
    Usage:
    input: xx_tab, ff_tab: two list of floats with the same length
        xx_tab: list of values, usually with a constant delta_x step size, 
                example: xx_tab=[0.0, 0.1, 0.2, 0.3, ...]
        ff_tab: values of a function f(x) at the values of xx_tab
                example: ff_tab=[f(0.0), f(0.1), f(0.2), ...]
        F0: float value
                initial value of integral
    output: intff : approximate integral of f(x) function given in the lists with F0 initial value.
        intff = [F0, F0+∫_0.0^0.1 f(x) dx,  F0+∫_0.0^0.2 f(x) dx, ....]

    Example usage: calculating position from velocity and original position data.
                x_tab=integral(t_tab, vx_tab, x0)
    '''
    N=len(xx_tab)     # size of the array
    intff=[0.0]*N     # array for integral values
    intff[0]=F0       # start value
    for i in range(1,N):
        intff[i]=intff[i-1]+(xx_tab[i]-xx_tab[i-1])*(ff_tab[i]+ff_tab[i-1])/2.0
    return(intff)

# Fill a list with a series of numbers

def fill_list_series(value_0, value_max, delta_value)->list:
    """     
    Usage:
    input: value_0, value_max, delta_value: floating point numbers
    output: series: list of floats
            series=[value_0, value_0+1*delta_value, value_0+2*delta_value, ...]

    Example usage: fill table of time values.
                t_tab=fill_list_series(0.0, 10.0, 0.1)
                output: t_tab=[0.0, 0.1, 0.2, .... , 9.9, 10.0]
    """
    series = [] # start with empty list
    value = value_0
    i = 0
    while (value <= value_max + delta_value/100):
        value = value_0 + i*delta_value
        series.append(value)
        i+=1
    
    return series

# Calculate function for every element of a list

def calc_func_list(value_list:list, func )->list:
    """ 
    Usage:
    input: value_list: list of floats
        func : a function name
    output: f_list: list of floats
            f_list=[func(value_list[0]), func(value_list[1]), func(value_list[2]), ... ]

    Example usage: fill a table for a function 
            x_tab=calc_func_list(t_tab, math.sin)
        output=[math.sin(t_tab[0]), math.sin(t_tab[1]), math.sin(t_tab[2]), ...]
    """    
    f_list=[] # start with empty list
    
    for value in value_list:   # loop over the list
        f_list.append(func(value))  # add a new element to the result list
    
    return f_list

# find index in a table (list), where a specific value is reached

def find_ind(tab:list, value)->list:
    """ 
    Usage:
    input: tab: list of floats
        value: value to search in tab
    output: ind_list: list of integers
            for every element 'ind' of ind_list we know, that value is between tab[ind] nad tab[ind+1]

    Examle usage: find the place in the position list, where the body reached a specific postition
            reached_10m=find_int(x_tab, 10.0)
        After this reached_10m will contain the index values where x_tab contains 10.0
        This list can be empty, contain 1 or more elements.
    """
    ind_list = []  # list of index values, where value is between tab[ind] and tab[ind+1]
    for ind in range(len(tab)-1):
        if (tab[ind] == value):
            ind_list.append(ind)
        elif (tab[ind] < value) and (tab[ind+1] > value):   # reach value from small values:
            ind_list.append(ind)
        elif (tab[ind] > value) and (tab[ind+1] < value):   # reach x_search from large values:
            ind_list.append(ind)
            
    return ind_list

# find index of maximum value in table

def find_max_ind(tab:list)->int:
    """ 
    Usage:
    input: tab: list of floats
    output: ind_max: integer
            tab[ind_max] is the (first) maximum of tab values

    Example usage: find maximum position during a motion.
            ind_top=find_max(x_tab)

    Note: This is the index, not the maximum value! x_top=x_tab[ind_top]
    """
    ind_max = 0
    value_max = tab[ind_max]
    
    for ind in range(1, len(tab)):
        value = tab[ind]
        if value > value_max:   # new maximum found
            value_max = value
            ind_max = ind
            
    return ind_max


# find index of minimum value in table

def find_min_ind(tab:list)->int:
    """ 
    Usage:
    input: tab: list of floats
    output: ind_min: integer
            tab[ind_min] is the (first) minimum of tab values

    Example usage: find minimum position during a motion.
            ind_bottom = find_min(x_tab)

    Note: This is the index, not the minimum value! x_bottom = x_tab[ind_bottom]
    """
    ind_min = 0
    value_min = tab[ind_min]
    
    for ind in range(1, len(tab)):
        value = tab[ind]
        if value < value_min:   # new minimum found
            value_min = value
            ind_min = ind
            
    return ind_min

## 1D functions
# calculate displacement from the beginning, assuming 1D motion

def calc_displacement_1D(x_tab:list)->list:
    """ Usage:
    input: x_tab
        table of x coordinates
    output: displ_tab 
            table of displacements relative to the start point

    Example usage: calculate displacement from the start point
            displ_tab=calc_displacement(x_tab)
    """
    N = len(x_tab)     # size of the array
    displ_tab = [0.0]*N  
    for i in range(1, N):
        displ_tab[i] = x_tab[i]-x_tab[0]
        
    return displ_tab
        
# calculate distance covered from the beginning, assuming 1D motion

def calc_distance_covered_1D(x_tab: list) -> list:
    """ Usage:
    input: x_tab
        table of x coordinates
    output: distcov_tab 
            table of total distance covered from the start point

    Example usage: calculate distance coevered from the start point
            dist_tab=calc_distance_coveredt(x_tab) """
    N = len(x_tab)     # size of the array
    distcov_tab = [0.0]*N  
    for i in range(1,N):
        distcov_tab[i] = distcov_tab[i-1]+abs(x_tab[i]-x_tab[i-1])
        
    return distcov_tab

## 2D functions
# calculate magnitude from two components

def calc_abs_2D(xc_tab: list, yc_tab: list)->list:
    '''
    Usage:
    input: xc_tab, yc_tab:
       x and y components of the vector
    Example usage: calculate speed from velocity components
       vabs_tab = calc_abs_2D(vx_tab, vy_tab)
    '''
    N = len(xc_tab)     # size of the array
    abs_tab = [0.0]*N 
    
    for i in range(N):
        abs_tab[i] = (xc_tab[i]**2 + yc_tab[i]**2)**0.5
        
    return abs_tab


# calculate displacement from the beginning, assuming 2D motion

def calc_displacement_2D(x_tab:list, y_tab:list):
    '''
    Usage:
    input: x_tab, y_tab
       tables of x and y coordinates
    output: displ_x_tab, displ_y_tab 
        table of displacements relative to the start point

    Example usage: calculate displacement from the start point
        displ_x_tab, displ_y_tab = calc_displacement_2D(x_tab, y_tab)
    '''
    displ_x_tab = calc_displacement_1D(x_tab)
    displ_y_tab = calc_displacement_1D(y_tab)
        
    return displ_x_tab, displ_y_tab

# calculate distance covered from the beginning, assuming 2D motion

def calc_distance_covered_2D(x_tab:list, y_tab:list)->list:
    '''
    Usage:
    input: x_tab, y_tab
       table of x and y coordinates
    output: distcov_tab 
        table of total distance covered from the start point

    Example usage: calculate distance coevered from the start point
        dist_tab = calc_distance_covered_2D(x_tab,y_tab)
    '''
    N = len(x_tab)     # size of the array
    distcov_tab = [0.0]*N  
    for i in range(1, N):
        delta_x = x_tab[i]-x_tab[i-1]
        delta_y = y_tab[i]-y_tab[i-1]
        distcov_tab[i] = distcov_tab[i-1] + (delta_x**2+delta_y**2)**0.5
        
    return distcov_tab

