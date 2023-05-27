class Vector(list):
    def __getitem__(self, i):
        if isinstance(i, slice):
            return Vector(list.__getitem__(self, i))
        
        return list.__getitem__(self, i)

    def __sub__(self, other):
        if not isinstance(other, Vector):
            raise TypeError(f"TypeError: unsupported operand type(s) for +: 'Vector' and '{type(other)}'")

        if len(self) != len(other):
            raise TypeError("TypeError: can't add vectors with different sizes")

        return Vector((s - other[i] for i, s in enumerate(self)))
    
    def __mul__(self, other):
        if not isinstance(other, Vector):
            raise TypeError(f"TypeError: unsupported operand type(s) for *: 'Vector' and '{type(other)}'")

        if len(self) != len(other):
            raise TypeError("TypeError: can't mul vectors with different sizes")
        
        return Vector((s * other[i] for i, s in enumerate(self)))

    def mul_by_coef(self, coef):
        if not (isinstance(coef, float) or isinstance(other, int)):
            raise TypeError(f"TypeError: unsupported operand type(s) for mul_by_coef: 'Vector' and '{type(coef)}'")

        return Vector((s * coef for s in self))
    
    def div_by_coef(self, coef):
        if not (isinstance(coef, float) or isinstance(coef, int)):
            raise TypeError(f"TypeError: unsupported operand type(s) for truediv_by_coef: 'Vector' and '{type(coef)}'")

        return Vector((s / coef for s in self))
        

if __name__ == '__main__':
    # read input data
    inp = input().strip()
    n, m = (int(i) for i in inp.split())

    matrix_ab = []
    for _ in range(n):
        inp = input().strip()
        matrix_ab.append(Vector(float(i) for i in inp.split()))

    # forward Gauss, transform matrix to echelon form
    k = 0
    for i, ab_i in enumerate(matrix_ab):
        if k == m:
            break

        if ab_i[k] == 0:
            continue
        
        ab_i = ab_i.div_by_coef(ab_i[k])
        
        for j in range(i + 1, n):
            matrix_ab[j] = matrix_ab[j] - ab_i.mul_by_coef(matrix_ab[j][k])
        
        k += 1
    
    # compute ranks
    rank_ab = sum(any(ab) for ab in matrix_ab)
    rank_a = sum(any(ab[:-1]) for ab in matrix_ab)
    if rank_ab != rank_a:
        print("NO")
    elif rank_ab < m:
        print("INF")
    else:
        print("YES")
        # reverse Gauss
        ans = Vector((1,) * m)
        k = m - 1
        for ab in reversed(matrix_ab):
            a = ab[:-1]
            b = ab[-1]

            a = a * ans
            b -= sum(a[k+1:])

            ans[k] = b / a[k]

            k -= 1
        print(' '.join(str(v) for v in ans))
