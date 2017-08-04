# complementary R code for "Extracting active modules from multilayer PPI network: a continuous optimization approach"
projectedgradient <- function(W,z,lambda=1,maxiter=100){
    k <- colSums(W)
    D <- diag(k)
    L <- D-W
    n = dim(W)[1]
    eta1 = rep(1,n)

    #x0 = rep(1/n,n)
    
    x0=runif(n)
    x0 <- sum(eta1)*x0-sum(x0)*eta1
    x0 <- x0/sqrt(sum(x0^2))*sqrt(n)

    x = x0
    epsilon = 1e-9
    beta = 0.1
    
    grad = L%*%x-lambda*z
    f_x = 0.5*(t(x)%*%L%*%x)[1,1]-lambda*(t(z)%*%x)
    func = numeric(maxiter)
    x_cand = x
    for (iteration in 1:maxiter){
        func[iteration] = f_x
        x_cand <- x-beta*grad    
        # orthogonalization: sum(\|x\|)=0, construct orthogonal vector
        x_cand <- sum(eta1)*x_cand-sum(x_cand)*eta1

        # normalization: \|x\|_2^2=n
        x_cand <- x_cand/sqrt(sum(x_cand^2))*sqrt(n)

        #print(sum(abs(x_cand-x)^2)^(1/2))
        if(sum(abs(x_cand-x)^2)^(1/2) < epsilon){
            break
        }
        x=x_cand
        grad = L%*%x-lambda*z
        f_x = 0.5*t(x)%*%L%*%x-lambda*(t(z)%*%x)
    }
    func <- func[1:iteration]
    return(list(func,x))
}

ModuleExtraction <-function(g, scores, maxsize){
    A <- as_adjacency_matrix(g)
    W <- as.matrix(A)
    z <- scores
    xx <- projectedgradient(W,z,lambda=1,maxiter=100)
    idxp <- which(xx[[2]]>0)
    subg <- subgraph(g, idxp)
    cp <- components(subg)
    grps <- groups(cp)
    largm <- which(cp$csize==max(cp$csize))[1]
    giancomp <- which(cp$membership==largm)
    idxp <- match(names(giancomp), rownames(W))
    while (length(giancomp) > maxsize) {
        W <- W[idxp,idxp]
        z <- z[idxp]
        tmpg <- graph_from_adjacency_matrix(W,mode='undirected')
        xx <- projectedgradient(W,z,lambda=1,maxiter=100)
        idxp <- which(xx[[2]]>0)
        subgnames <- rownames(W)[idxp]
        idxp <- match(subgnames,names(z))
        subg <- subgraph(tmpg, idxp)
        cp <- components(subg)
        largm <- which(cp$csize==max(cp$csize))[1]
        giancomp <- which(cp$membership==largm)
        idxp <- match(names(giancomp), rownames(W))
    }
    return(subg)
}

pgLineSearch <- function(W,z,lambda=1,maxiter=100){
    k <- colSums(W)
    D <- diag(k)
    L <- D-W
    n = dim(W)[1]
    x = x0
    epsilon = 1e-9
    beta = 0.1
    eta1 = rep(1,n)
    grad = L%*%x-lambda*z
    f_x = 0.5*(t(x)%*%L%*%x)[1,1]-lambda*(t(z)%*%x)
    func = numeric(maxiter)
    x_cand = x
    for (iteration in 1:maxiter){
        func[iteration] = f_x
        alpha <- LineSearch(x,grad)
        print(alpha)
        x_cand <- x-alpha*grad    
        # orthogonalization: sum(\|x\|)=0, construct orthogonal vector
        x_cand <- sum(eta1)*x_cand-sum(x_cand)*eta1

        # normalization: \|x\|_2^2=n
        x_cand <- x_cand/sqrt(sum(x_cand^2))*sqrt(n)

        #print(sum(abs(x_cand-x)^2)^(1/2))
        if(sum(abs(x_cand-x)^2)^(1/2) < epsilon){break}
        x=x_cand
        grad = L%*%x-lambda*z
        f_x = 0.5*(t(x)%*%L%*%x)[1,1]-lambda*(t(z)%*%x)
    }
    func <- func[1:iteration]
}

# search step size alpha
LineSearch <- function(x,grad){
    n = length(x)
    alpha = 1
    beta = 0.1
    sigma = 0.01
    eta1 = rep(1,n)
    for (inner_iter in 1:20) {
        xp = x-alpha*grad
        xp = sum(eta1)*xp-sum(xp)*eta1
        xp = xp/sqrt(sum(xp^2))*sqrt(n)
        
        dx = xp-x
        suff_decr = 0.5*(t(dx)%*%L%*%dx)[1,1]/n-lambda*(t(z)%*%dx)-sigma*t(grad)%*%dx < 0

            if (inner_iter==1)
                decr_alpha = !suff_decr
            if(decr_alpha){
                if (suff_decr){
                    break
                } else {
                    alpha = alpha * beta
                }
            } else {
                if (!suff_decr) {
                    break
                } else{
                    alpha = alpha/beta
                }
            }
    }
    return (alpha)
}
