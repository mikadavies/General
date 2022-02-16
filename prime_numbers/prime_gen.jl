function generator_1(max=10000)
    primes = [2, 3]
    function test(n)
        d = 2
        while d <= n/2
            if n%d != 0
                d += 1
            else
                return false
            end
        end
        if d > n/2
            return true
        end
    end
    for n in 4:max
        if test(n) == true
            append!(primes, n)
        end
    end
    return primes
end

function generator_2(max=10000)
    p = [2, 3]
    function test(n)
        d = 2
        i = 1
        while d <= n/2
            if n%d != 0
                if d < last(p)
                    d = p[i]
                    i += 1
                else
                    d += 1
                end
            else
                return false
            end
        end
        if d > n/2
            return true
        end
    end
    for n in 4:max
        if test(n) == true
            append!(p, n)
        end
    end
    return p
end

function generator_3(max=10000)
    p = [2, 3]
    function test(n)
        i = 1
        while i <= length(p)
            d = p[i]
            if n%d != 0
                i += 1
            else
                return false
            end
        end
        if i > length(p)
            return true
        end
    end
    for n in 4:max
        if test(n) == true
            append!(p, n)
        end
    end
    return p
end

function generator_4(max=10000)
    p = [2, 3]
    for n in 4:max
        i = 1
        while i <= length(p)
            if n%p[i] != 0
                if i == length(p)
                    append!(p, n)
                else
                    i+=1
                end
            else
                break
            end
        end
    end
    return p
end

function generator_5(max=10000)
    p = [2, 3]
    function test(n)
        d = 2
        i = 1
        while d <= n/2
            if n%d != 0
                if d <= last(p) #|| d == last(p)
                    d = p[i]
                    i += 1
                #else break
                end
            else return false
            end
        end
    end
    for n in 4:max
        if test(n) != false
            append!(p, n)
        end
    end
    return p
end

function generator_6(max=10000)
    p = [2, 3]
    n = 5
    while n < max
        i = 1
        while i <= length(p)
            if p[i] > n/2
                append!(p, n)
                break
            elseif n%p[i] != 0
                i+=1
            else
                break
            end
        end
        n += 2
    end
    return p
end

function generator_7(max=10000) # AKA Erastothenes
    p = [x for x in 2:max if isodd(x)]
    for n in p
        i = n
        while i <= max/n
            m = i*n
            deleteat!(p, findall(x->x==m, p))
            i += 1
        end
    end
    return append!([2] ,p)
end

function generator_8(max=10000) # AKA Atkin?
    list = [false for n in 1:max]
    p = [2,3,5]
    function alg1()
        r1 = [1, 13, 17, 29, 37, 41, 49, 53]
        x = 1
        y = 1
        while x <= max
            while y <= max
                n = 4*x^2+y^2
                if n%60 in r1 && n <= max
                    list[n] = true
                end
                y += 2
            end
            x += 1
        end
    end
    function alg2()
        r2 = [7, 19, 31, 43]
        x = 1
        y = 2
        while x <= max
            while y <= max
                n = 3*x^2+y^2
                if n%60 in r2 && n <= max
                    list[n] = true
                end
                y += 2
            end
            x += 2
        end
    end
    function alg3()
        r3 = [11, 23, 47, 59]
        x = 2
        N = 1
        while x <= max
            while N < x
                y = x-N
                n = 3*x^2-y^2
                if n%60 in r3 && n <= max
                    list[n] = true
                end
                N += 2
            end
            x += 1
        end
    end
    alg1()
    alg2()
    alg3()
    n = 7
    while n != nothing && n < length(list)
        n = findnext(list, n)
        if n != nothing
            append!(p, n)
            s = n^2
            n += 1
            if s <= max
                for i in 1:Int(floor(max/s))
                    np = s*i
                    list[np] = false
                end
            end
        end
    end
    return p
end
