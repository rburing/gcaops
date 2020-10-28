def selection_sort(lst):
    sign = 1
    for i in range(0, len(lst)-1):
        j = min(range(i, len(lst)), key=lst.__getitem__)
        if i != j:
            lst[i],lst[j] = lst[j],lst[i]
            sign *= -1
    return sign
