def selection_sort(lst):
    sign = 1
    for i in range(0, len(lst)-1):
        j = min(range(i, len(lst)), key=lst.__getitem__)
        if i != j:
            lst[i],lst[j] = lst[j],lst[i]
            sign *= -1
    return sign

def selection_sort_graded(lst, degrees):
    sign = 1
    for i in range(0, len(lst)-1):
        j = min(range(i, len(lst)), key=lst.__getitem__)
        if i != j:
            if (degrees[i]*sum(degrees[k] for k in range(i+1,j+1)) + degrees[j]*sum(degrees[k] for k in range(i+1,j))) % 2 == 1:
                sign *= -1
            lst[i],lst[j] = lst[j],lst[i]
            degrees[i],degrees[j] = degrees[j],degrees[i]
    return sign
