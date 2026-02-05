###Exemplo 5.11 do Incropera - Condução transiente. 

#Condução transiente em parede plana com 2L = 20 mm. 
#Balanço de energia discretizado no interior da placa plana (condução), adiabática na parede interna e convecção na parede externa. 

using Plots
using DataFrames

m = 5
n = 2000 #Quantidade de passos de tempo
p = 0:1:n
T_array = [Vector{Float64}() for i in 0:n]
#Propriedades dadas: 
L = 10e-3
h = 1100
T_inf = 250 + 273.15 
q_inicial = 1e7
q_final = 2e7

t_final = 1.5 
k = 30 
α = 5e-6

#Parâmetros: 

Δx = L/5 #Passo entre os nós.  
Bi = h*Δx/k #Número de Biot. 

#Critério de estabilidade: Fo(1 + Bi) <= 1/2

Fo = 1/2 * 1/(1 + Bi) 
println("Fo = ", Fo)

#Passo de tempo: 

Δt_estabilidade = Fo*Δx^2/α #Fo = α*t/L^2 

println("Passo de tempo: ", Δt_estabilidade, " s")
#Escolher Δt <= Δt_estabilidade

Δt = 0.3 #Critério de estabilidade (Δt <= Δt_estabilidade)

Fo_new = α * Δt / (Δx)^2

#Valores iniciais para distribuição de temperatura ao longo da placa plana: T(x) = 16.67*(1 - x^2/L^2) + 340.91 [°C]

L_list = []
t_list = [p*Δt for p in collect(p)]
for i in collect(0:m)
    push!(L_list, L*i/m)
    push!(T_array[1], 16.67 * (1 - L_list[i+1]^2/L^2) + 340.91)
end

for i in 0:length(collect(p))-2, j in collect(0:m)
    #j => nó, i => discretização do tempo (número inteiro)
    if j == 0
        push!(T_array[i+2], 0.375*(2*T_array[i+1][2] + 2.67) + 0.250 * T_array[i+1][1])
    elseif j == 5
        push!(T_array[i+2], 0.750*(T_array[i+1][5] + 19.67) + 0.195*T_array[i+1][6])
    else
        push!(T_array[i+2], 0.375*(T_array[i+1][j] + T_array[i+1][j+2] + 2.67) + 0.250*T_array[i+1][j+1])
    end
end

df = DataFrame(reduce(vcat, permutedims.(T_array)), [:T₀, :T₁, :T₂, :T₃, :T₄, :T₅])
insertcols!(df, 1, :p => collect(p))
insertcols!(df, 2, :Δt => t_list)

#Permutedims => transforma vector (que é uma columa) em uma matrix transposta do vector. 
#vcat => concatena verticalmente a nova matrix, se usado o splatter (...).
#reduce => apliica uma função (vcat) em um conjunto de dados.

T0_list = [T_array[i][1] for i in 1:length(T_array)]
T5_list = [T_array[i][end] for i in 1:length(T_array)]

ErrorT0 = minimum(filter(!isnothing, [abs(T0_list[i+1] - T0_list[i]) <= 1e-4 ? T0_list[i] : nothing for i in 1:length(T0_list)-1]))
ErrorT5 = minimum(filter(!isnothing, [abs(T5_list[i+1] - T5_list[i]) <= 1e-4 ? T5_list[i] : nothing for i in 1:length(T5_list)-1]))


T_plot = plot(t_list, T0_list, xlabel = "t, s", ylabel = "T, °C", label = "Inner wall temperature", linestyle = :dot, color = :black)
plot!(T_plot, t_list, T5_list, label = "Outer wall temperature", color = :black)
plot!(t_list, [ErrorT0 for i in 1:length(t_list)], linestyle = :dash, color = :grey, label = "")
plot!(t_list, [ErrorT5 for i in 1:length(t_list)], linestyle = :dash, color = :grey, label = "")

ticks = [320, 360, 400, 440, 480]
targets = [round(ErrorT0[1],digits = 2), round(ErrorT5[1], digits = 2)]
all_postions = sort([ticks..., targets...])
#Splatter para quebrar a lista em argumentos, variáveis. 
yticks!(T_plot, all_postions)
2savefig(T_plot, "C:/Users/gabri/Desktop/Incropera511.png")