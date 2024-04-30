\name{POETKhat}
\alias{POETKhat}
\title{Khat - number of factors in approximate factor model

}
\description{This function is for calculating the optimal number of factors
in an approximate factor model.
}

\usage{
POETKhat(Y)
}

\arguments{
  \item{Y}{
p by n matrix of raw data, where p is the dimensionality, n
is the sample size. It is recommended that Y is de-meaned,
i.e., each row has zero mean.
}
}

\details{ This method was proposed by Bai & Ng (2002) and Hallin & Liska (2007).
They propose two penalty functions and in turn minimize the corresponding
information criteria. Notice that this method may underestimate K. POET is very robust to over-estimating K. But under-estimating K
can result to VERY BAD performance. Therefore we strongly recommend
choosing a relatively large K (normally less than 8) to avoid
missing any important common factor.


}

\value{

\item{K1HL}{estimated number of factors based on the first infomation criterion
using Hallin & Liska method}
 \item{K2HL}{estimated number of factors based on the second information
criterion using Hallin & Liska method}
\item{K1BN}{estimated number of factors based on the first infomation criterion
using Bai & Ng method}
 \item{K2BN}{estimated number of factors based on the second information
criterion using Bai & Ng method}
}

\references{
Bai,Ng,2002.Determining the number of factors in approximate factor models.
Econometrica 70,191-221.

Hallin,Liska,2007.Determining the number of factors in the general dynamic 
factor model.JASA 102,603-617.

Alessi,Barigozzi,Capasso,2010. Improved penalization for determining the number
of factors in approximate factor models. Statistics and Probability Letters 80,
1806-1813.
}


\author{Jianqing Fan, Yuan Liao, Martina Mincheva

}

\examples{
p=100
n=100
Y<-array(rnorm(p*n),dim=c(p,n))
K<-POETKhat(Y)
}

