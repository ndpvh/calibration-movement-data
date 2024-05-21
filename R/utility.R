# A vectorized sequence function
multi_seq <- Vectorize(seq.default, 
                       vectorize.args = c("from", "to", "by", "length.out"))