'''
tasks:
restrict to binary x-rays.
find holes. check possible arc length combinations/x-rays and try to form a permutation.
'''
from a import get_x_ray_strings_for as albert


def generate_binary_x_rays(n):
    tableu_uses = [0]
    def get_ways(length, num_ones, weighted_sum, tableu, depth=0):
        '''
        for a given length of binary number, number of ones, and weighted sum (starting at index 0 on the left),
        return a list of satisfying binary numbers. 
        Use the tableu to memoize results.
        e.g. tableu[(3, 2, 2)] = ['101'] (length 3, 2 ones, sum = 0 + 2 = 2 and this is the only way).
        '''
        if length < 0 or num_ones < 0 or weighted_sum < 0:
            return [] # no ways to form such an x-ray
        if length == 0 and num_ones == 0 and weighted_sum == 0:
            return ['']
        if (length, num_ones, weighted_sum) in tableu:
            tableu_uses[0] += 1
            return tableu[(length, num_ones, weighted_sum)]
        if (num_ones * (length - num_ones) < weighted_sum):
            # the maximum backwards weighted sum of the remaining x-ray is not enough.
            # we dont have enough space left to possibly make the weighted sum at each index at least the number of ones squared. 
            return []

        # if we start with a 1, the sum of the remaining x-ray is 1 more to balance it out (but all indices decrease, so subtract num_ones)
        ways_start_1 = get_ways(length-1, num_ones-1, weighted_sum-num_ones+1, tableu, depth+1)
        # if we start with a 0, the weighted_sum of the remaining x-ray is the same (but all indices decrease, so subtract num_ones)
        ways_start_0 = get_ways(length-1, num_ones, weighted_sum-num_ones, tableu, depth+1)

        ways_start_0 = [f'0{x}' for x in ways_start_0]
        ways_start_1 = [f'1{x}' for x in ways_start_1]

        total_ways = ways_start_1 + ways_start_0
        tableu[(length, num_ones, weighted_sum)] = total_ways

        if depth == 0:
            print("Number of entries saved in tableu", len(tableu.keys()))
            print("Total tableu x-rays: ", sum([len(ways) for ways in tableu.values()]))
            print("Tableu uses: ", tableu_uses[0])
            print("Number of x-rays: ", len(total_ways))
        return total_ways
        
    tableu = {}
    return get_ways(2*n-1, n, n*(n-1), tableu)
        

# def count_binary_x_rays(n): 
#     '''
#     Don't actually generate x-rays, but count the number of x-rays for a given n.
#     '''
#     def get_ways(length, num_ones, sum, tableu, depth=0):
#         if length < 0 or num_ones < 0 or sum < 0:
#             return 0
#         if length == 0:
#             return 1 if num_ones == 0 and sum == 0 else 0
#         if (length, num_ones, sum) in tableu:
#             return tableu[(length, num_ones, sum)]
#         else:
#             # recursive step: vc
#             ways_start_1 = get_ways(length-1, num_ones-1, sum-num_ones+1, tableu, depth+1)
#             ways_start_0 = get_ways(length-1, num_ones, sum-num_ones, tableu, depth+1)
#             tableu[(length, num_ones, sum)] = ways_start_1 + ways_start_0
#             if depth == 0:
#                 print("Size of counting tableu: ", len(tableu.keys()))
#             return ways_start_1 + ways_start_0
        
#     tableu = {(0, 0, 0): 1} # map tuples of the form (length, num_ones, sum) to number of ways to form such an x-ray
#     return get_ways(2*n-1, n, n*(n-1), tableu)



for n in range(6, 7):
    print("\nn = ", n)

    x_rays = set(generate_binary_x_rays(n))

    for x_ray in x_rays:
        print(x_ray)


# 1100011 -> (1, 2, 6, 7) sum = 16 = 4^2
# 1 + 2 = 3 not >= 2^2...



'''
110001101
101010101
101001110
100111001
100110110
011100101
011011001
011010110
010111010
001111100

101010101
101001110
100111001
100110110
011100101
011011001
011010110
010111010
001111100

'''