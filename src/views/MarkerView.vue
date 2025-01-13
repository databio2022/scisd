<template>
  <div class="container mx-auto mt-4">
    <div class="grid grid-cols-12 gap-4">
      <div class="col-span-12 p-4">
        <div class="card bg-base-100 border">
          <div class="p-4 bg-gray-200">
            <h2 class="text-lg font-bold">Cell Clustering And Gene Expression</h2>
          </div>
          <div class="card-body">
            <div class="flex justify-center space-x-10">
              <div class="relative">
                <label class="input input-bordered flex items-center gap-2 max-w-sm">
                  <input
                    type="text"
                    class="grow"
                    placeholder="Search gene, e.g. ACTB"
                    v-model="selectedGene"
                  />
                  <svg
                    xmlns="http://www.w3.org/2000/svg"
                    viewBox="0 0 16 16"
                    fill="currentColor"
                    class="h-4 w-4 opacity-70"
                  >
                    <path
                      fill-rule="evenodd"
                      d="M9.965 11.026a5 5 0 1 1 1.06-1.06l2.755 2.754a.75.75 0 1 1-1.06 1.06l-2.755-2.754ZM10.5 7a3.5 3.5 0 1 1-7 0 3.5 3.5 0 0 1 7 0Z"
                      clip-rule="evenodd"
                    />
                  </svg>
                </label>
                <div
                  v-if="!isGeneSelected && filteredGenes.length"
                  class="absolute top-full left-0 w-full z-10 bg-white border shadow-md"
                >
                  <VirtualList :list="filteredGenes" @gene-selected="handleGeneSelected" />
                </div>
              </div>
              <button class="btn btn-accent" @click="searchgene">Search</button>
            </div>
            <div class="flex justify-center">
              <!-- waiting animation-->
              <div class="flex items-center justify-center m-24 p-24" v-if="isLoading">
                <div
                  class="animate-spin rounded-full h-12 w-12 border-t-4 border-blue-500 border-opacity-75"
                ></div>
              </div>

              <div v-if="loading">
                <div class="flex justify-between">
                  <img :src="imgurl" class="w-1/2 h-[400px]" />
                  <BoxplotChart :data="boxPlotData" :gene="selectedGene" />
                </div>
              </div>
            </div>
          </div>
        </div>
        <div class="card bg-base-100 border mt-4">
          <div class="p-4 bg-gray-200">
            <h2 class="text-lg font-bold">Differentially Expressed Genes</h2>
          </div>
          <div class="card-body">
            <div>
              <!-- waiting animation-->
              <div class="flex items-center justify-center m-24 p-24" v-if="isLoading">
                <div
                  class="animate-spin rounded-full h-12 w-12 border-t-4 border-blue-500 border-opacity-75"
                ></div>
              </div>

              <div v-if="loading">
                <div class="flex space-x-10">
                  <div>
                    <h1 class="font-bold py-4">Select Cell Type:</h1>
                    <select class="select select-bordered w-full max-w-xs mb-4" v-model="celltype">
                      <option v-for="type in celltypes" :key="type" :value="type">
                        {{ type }}
                      </option>
                    </select>
                  </div>
                  <div>
                    <h1 class="font-bold py-4">Select Comparison Groups:</h1>
                    <select
                      class="select select-bordered w-full max-w-xs mb-4"
                      v-model="comparison"
                    >
                      <option value="cell">Cell type-specific genes</option>
                      <option v-for="type in multigroups" :key="type" :value="type">
                        {{ type }}
                      </option>
                    </select>
                  </div>
                  <div>
                    <h1 class="font-bold py-4 text-white">aaa</h1>
                    <button class="btn btn-accent" @click="getdeg">Search</button>
                  </div>
                </div>
                <div v-if="showdeg">
                  <TheDegDataTable :data="degtabledata"></TheDegDataTable>
                </div>
                <div v-else><TheDataTable :data="tabledata"></TheDataTable></div>
                <hr />
                <div>
                  <div>
                    <h1 class="font-bold py-4">
                      GO and KEGG Analysis<span class="text-red-500"
                        >(Please click 'Search' button to update data when above table
                        changes)</span
                      >:
                    </h1>
                    <div class="flex justify-between">
                      <div class="flex space-x-4">
                        <select class="select select-bordered w-full max-w-xs mb-4" v-model="go">
                          <option value="KEGG">KEGG</option>
                          <option value="BP">Biological Process</option>
                          <option value="CC">Cellular Component</option>
                          <option value="MF">Molecular Function</option>
                        </select>
                        <button class="btn btn-accent" @click="getgo">Search</button>
                      </div>
                      <!-- Open the modal using ID.showModal() method -->
                      <button class="btn btn-accent" @click="showModalAndPlot">Top10 Plot</button>
                      <dialog id="my_modal_1" class="modal" ref="myModal">
                        <div class="modal-box max-w-4xl">
                          <h3 class="text-lg font-bold">
                            The Figure show top 10 significant terms
                          </h3>
                          <div v-if="showbar">
                            <BarPlotComponent :data="xbardata" :layout="xlayout"></BarPlotComponent>
                          </div>
                          <div class="modal-action">
                            <form method="dialog">
                              <!-- if there is a button in form, it will close the modal -->
                              <button class="btn" @click="fclose">Close</button>
                            </form>
                          </div>
                        </div>
                      </dialog>
                    </div>
                    <TheGoDataTable :data="tabledata1"></TheGoDataTable>
                  </div>
                </div>
              </div>
            </div>
          </div>
        </div>
      </div>
    </div>
  </div>
  <p>{{ id }}</p>
</template>

<script>
import TheDataTable from '../components/TheDataTable.vue'
import TheGoDataTable from '../components/TheGoDataTable.vue'
import BoxplotChart from '../components/BoxplotChart.vue'
import axios from 'axios'
import VirtualList from '../components/VirtualList.vue'

import TheDegDataTable from '@/components/TheDegDataTable.vue'
import BarPlotComponent from '../components/BarplotChart.vue'

export default {
  components: {
    TheDataTable,
    TheGoDataTable,
    BoxplotChart,
    VirtualList,
    TheDegDataTable,
    BarPlotComponent
  },
  data() {
    return {
      loading: false,
      isLoading: true,
      imgurl: null,

      groups: [],

      selectedGSE: this.$route.params.id,
      selectedGene: 'ACTB',

      showdeg: false,
      degtabledata: [],
      tabledata: [],
      tabledata1: [],
      boxPlotData: [],

      umapData: [],
      substatus: false,
      celltype: '',
      celltypes: [],
      comparison: 'cell',
      multigroups: [],
      go: 'KEGG',
      genes: [],
      isGeneSelected: false,
      xbardata: [],
      xlayout: null,
      showbar: false
    }
  },
  mounted() {
    this.fetchData()
  },
  computed: {
    filteredGenes() {
      if (!this.selectedGene || this.selectedGene === 'ACTB') {
        // 或者 if (this.selectedGene.trim() === '')
        return [] // 输入框为空时返回空数组
      }
      return this.genes.filter((gene) =>
        gene.toLowerCase().includes(this.selectedGene.toLowerCase())
      )
    }
  },

  watch: {
    selectedGene(newValue) {
      if (newValue === '') {
        this.isGeneSelected = false // 输入框为空时显示列表
      }
    }
  },

  methods: {
    handleGeneSelected(gene) {
      this.selectedGene = gene
      this.isGeneSelected = true
    },
    showModalAndPlot() {
      this.showbar = true
      this.$refs.myModal.showModal()
    },
    fclose() {
      this.showbar = false
    },

    searchgene() {
      this.isLoading = true
      this.loading = false
      axios
        .post('http://scisd.databio1.com/api/search_meta', {
          subject: this.selectedGSE,
          gene: this.selectedGene
        })
        .then((response) => {
          this.substatus = false
          // 处理服务器返回的数据
          //  console.log(JSON.parse(response.data[0].gex))
          //  console.log(JSON.parse(response.data[0].metadata))
          //  console.log(JSON.parse(response.data[0].umap_df))
          //  console.log(response.data[0])

          //  let data0 = JSON.parse(response.data[0].metadata)

          let data2 = JSON.parse(response.data[0].gex)

          // 在这里检查 data1.value 和 data2.value 是否已加载完成

          //  console.log(mergedData)

          //this.umapData = mergedData

          this.imgurl = 'http://scisd.databio1.com/api/images/' + this.selectedGSE

          //  console.log(this.umapData)
          // this.tabledata = JSON.parse(response.data[0].cor)
          //console.log(this.tabledata)

          //box
          let data1 = JSON.parse(response.data[0].metadata)
          const boxdata = data1.map((item, index) => ({
            ...item,
            x: data2[index].x
          }))

          // 获取所有 group 和 cluster
          const allClusters = [...new Set(boxdata.map((item) => item.celltype))]

          // 创建 Plotly 图表数据
          const newboxPlotData = allClusters.map((group) => {
            // 创建每个 group 的 trace 对象
            const trace = {
              y: [],
              x: [],
              name: group,
              type: 'box'
            }

            // 填充 y 和 x 数据

            const values = boxdata.filter((item) => item.celltype === group).map((item) => item.x)

            trace.y.push(...values) // 将该 cluster 的所有值添加到 y 数组
            trace.x.push(...Array(values.length).fill(group)) // 将对应数量的 cluster 添加到 x 数组

            return trace
          })

          //  console.log(newboxPlotData)
          this.boxPlotData = newboxPlotData
          this.isLoading = false
          this.loading = true
        })
    },
    getgo() {
      console.log(this.go)
      if (this.comparison.includes('vs')) {
        axios
          .get(
            'http://scisd.databio1.com/api/filterde_enrichment_' +
              this.go +
              '_table/' +
              this.selectedGSE +
              '/' +
              this.celltype.replace(/\s+/g, '_') +
              '/' +
              this.comparison
          )
          .then((response2) => {
            console.log('111')
            console.log(JSON.parse(response2.data))
            this.tabledata1 = JSON.parse(response2.data)
            const outputData = this.tabledata1
              .sort((a, b) => a['p.adjust'] - b['p.adjust']) // 按照 p.adjust 升序排序
              .slice(0, 10) // 取前10个元素
              .reduce(
                (acc, curr, index) => {
                  acc['x'].push(-Math.log10(curr['p.adjust']))

                  acc['y'].push(index)
                  acc['descriptions'].push(curr.Description)
                  acc['marker']['color'].push('rgba(54, 162, 235, 0.7)')
                  return acc
                },
                { x: [], y: [], marker: { color: [] }, orientation: 'h', descriptions: [] }
              )
            outputData['type'] = 'bar'

            this.xbardata = [outputData]
            console.log(this.xbardata)

            const layout = {
              yaxis: {
                autorange: 'reversed',
                automargin: true,
                tickmode: 'array',
                tickvals: outputData.y,
                ticktext: outputData.descriptions
              },
              xaxis: {
                title: '-log10(p.adjust)',
                automargin: true
              }
            }
            this.xlayout = layout
            this.showbar = true
          })
          .catch((error) => {
            console.error(error)
          })
      } else {
        axios
          .get(
            'http://scisd.databio1.com/api/cellmarker_enrichment_' +
              this.go +
              '_table/' +
              this.selectedGSE +
              '/' +
              this.celltype.replace(/\s+/g, '_')
          )
          .then((response2) => {
            console.log(JSON.parse(response2.data))
            this.tabledata1 = JSON.parse(response2.data)
            const outputData = this.tabledata1
              .sort((a, b) => a['p.adjust'] - b['p.adjust']) // 按照 p.adjust 升序排序
              .slice(0, 10) // 取前10个元素
              .reduce(
                (acc, curr, index) => {
                  acc['x'].push(-Math.log10(curr['p.adjust']))

                  acc['y'].push(index)
                  acc['descriptions'].push(curr.Description)
                  acc['marker']['color'].push('rgba(54, 162, 235, 0.7)')
                  return acc
                },
                { x: [], y: [], marker: { color: [] }, orientation: 'h', descriptions: [] }
              )
            outputData['type'] = 'bar'

            this.xbardata = [outputData]
            console.log(this.xbardata)

            const layout = {
              yaxis: {
                autorange: 'reversed',
                automargin: true,
                tickmode: 'array',
                tickvals: outputData.y,
                ticktext: outputData.descriptions
              },
              xaxis: {
                title: '-log10(p.adjust)',
                automargin: true
              }
            }
            this.xlayout = layout
            this.showbar = true
          })
          .catch((error) => {
            console.error(error)
          })
      }
    },
    getdeg() {
      console.log(this.comparison)
      if (this.comparison.includes('vs')) {
        axios
          .post('http://scisd.databio1.com/api/query_diffgene_group', {
            subject: this.selectedGSE,
            cluster: this.celltype.replace(/\s+/g, '_'),
            group: this.comparison
          })
          .then((response1) => {
            //  console.log(response1.data)
            this.degtabledata = response1.data
            this.showdeg = true
          })
          .catch((error) => {
            console.error(error)
          })
      } else {
        axios
          .post('http://scisd.databio1.com/api/query_diffgene', {
            subject: this.selectedGSE,
            cluster: this.celltype.replace(/\s+/g, '_')
          })
          .then((response1) => {
            //  console.log(response1.data)
            this.tabledata = response1.data
            this.showdeg = false
          })
          .catch((error) => {
            console.error(error)
          })
      }
    },

    async fetchData() {
      if (this.selectedGSE && this.selectedGene) {
        this.substatus = true
        this.loading = false
        this.isLoading = true
        try {
          const [response, responseb] = await Promise.all([
            axios.post('http://scisd.databio1.com/api/search_meta', {
              subject: this.selectedGSE,
              gene: this.selectedGene
            }),
            axios.get('http://scisd.databio1.com/api/gene_list/' + this.selectedGSE)
          ])

          this.substatus = false
          // 处理服务器返回的数据
          //  console.log(JSON.parse(response.data[0].gex))
          //  console.log(JSON.parse(response.data[0].metadata))
          //  console.log(JSON.parse(response.data[0].umap_df))
          // console.log(response.data)

          let data2 = JSON.parse(response.data[0].gex)

          //box
          let data1 = JSON.parse(response.data[0].metadata)
          const boxdata = data1.map((item, index) => ({
            ...item,
            x: data2[index].x
          }))

          // 获取所有 group 和 cluster
          const allClusters = [...new Set(boxdata.map((item) => item.celltype))]

          this.celltypes = allClusters
          this.celltype = allClusters[0]

          // 创建 Plotly 图表数据
          const newboxPlotData = allClusters.map((group) => {
            // 创建每个 group 的 trace 对象
            const trace = {
              y: [],
              x: [],
              name: group,
              type: 'box'
            }

            // 填充 y 和 x 数据

            const values = boxdata.filter((item) => item.celltype === group).map((item) => item.x)

            trace.y.push(...values) // 将该 cluster 的所有值添加到 y 数组
            trace.x.push(...Array(values.length).fill(group)) // 将对应数量的 cluster 添加到 x 数组

            return trace
          })

          //  console.log(newboxPlotData)
          this.boxPlotData = newboxPlotData
          this.imgurl = 'http://scisd.databio1.com/api/images/' + this.selectedGSE

          this.loading = true
          const [response1, response2, response3, response4] = await Promise.all([
            axios.post('http://scisd.databio1.com/api/query_diffgene', {
              subject: this.selectedGSE,
              cluster: this.celltype.replace(/\s+/g, '_')
            }),
            axios.get(
              'http://scisd.databio1.com/api/cellmarker_enrichment_KEGG_table/' +
                this.selectedGSE +
                '/' +
                this.celltype.replace(/\s+/g, '_')
            ),
            axios.get('http://scisd.databio1.com/api/vs_table/' + this.selectedGSE)
          ])

          this.tabledata = response1.data

          //  console.log(JSON.parse(response2.data))
          this.tabledata1 = JSON.parse(response2.data)
          const outputData = this.tabledata1
            .sort((a, b) => a['p.adjust'] - b['p.adjust']) // 按照 p.adjust 升序排序
            .slice(0, 10) // 取前10个元素
            .reduce(
              (acc, curr, index) => {
                acc['x'].push(-Math.log10(curr['p.adjust']))

                acc['y'].push(index)
                acc['descriptions'].push(curr.Description)
                acc['marker']['color'].push('rgba(54, 162, 235, 0.7)')
                return acc
              },
              { x: [], y: [], marker: { color: [] }, orientation: 'h', descriptions: [] }
            )
          outputData['type'] = 'bar'

          this.xbardata = [outputData]
          console.log(this.xbardata)

          const layout = {
            yaxis: {
              autorange: 'reversed',
              automargin: true,
              tickmode: 'array',
              tickvals: outputData.y,
              ticktext: outputData.descriptions
            },
            xaxis: {
              title: '-log10(p.adjust)',
              automargin: true
            }
          }
          this.xlayout = layout
          this.showbar = true
          const x = JSON.parse(responseb.data)
          const geneArray1 = x.map((item) => Object.values(item)[0])
          this.genes = geneArray1
          const x1 = JSON.parse(response3.data)
          const geneArray2 = x1.map((item) => Object.values(item)[1])
          //  this.genes = geneArray1
          this.multigroups = geneArray2
          //  console.log(x1)
        } catch (error) {
          console.error(error)
        } finally {
          this.isLoading = false
        }
      }
    }

    //    fetchData() {
    //      if (this.selectedGSE && this.selectedGene) {
    //        // 使用 axios 发送 POST 请求
    //        axios
    //          .post('http://scisd.databio1.com/api/search_meta', {
    //            subject: this.selectedGSE,
    //            gene: this.selectedGene
    //          })
    //          .then((response) => {
    //            this.substatus = false
    //            // 处理服务器返回的数据
    //            console.log(JSON.parse(response.data[0].gex))
    //            console.log(JSON.parse(response.data[0].metadata))
    //            console.log(JSON.parse(response.data[0].umap_df))
    //            console.log(response.data[0])
    //
    //          //  let data0 = JSON.parse(response.data[0].metadata)
    //          //  let data1 = JSON.parse(response.data[0].umap_df)
    //            let data2 = JSON.parse(response.data[0].gex)
    //
    //            // 在这里检查 data1.value 和 data2.value 是否已加载完成
    //          //  const mergedData = data1.map((item, index) => ({
    //          //    ...item,
    //          //    x: data2[index].x
    //          //  }))
    //
    //            //  console.log(mergedData)
    //
    //            //this.umapData = mergedData
    //
    //            this.imgurl = 'http://scisd.databio1.com/api/images/' + this.selectedGSE
    //
    //            //  console.log(this.umapData)
    //            // this.tabledata = JSON.parse(response.data[0].cor)
    //            //console.log(this.tabledata)
    //
    //            //box
    //            let data1 = JSON.parse(response.data[0].metadata)
    //            const boxdata = data1.map((item, index) => ({
    //              ...item,
    //              x: data2[index].x
    //            }))
    //
    //            // 获取所有 group 和 cluster
    //            const allClusters = [...new Set(boxdata.map((item) => item.celltype))]
    //
    //            this.celltypes = allClusters
    //            this.celltype = allClusters[0]
    //
    //            // 创建 Plotly 图表数据
    //            const newboxPlotData = allClusters.map((group) => {
    //              // 创建每个 group 的 trace 对象
    //              const trace = {
    //                y: [],
    //                x: [],
    //                name: group,
    //                type: 'box'
    //              }
    //
    //              // 填充 y 和 x 数据
    //
    //              const values = boxdata.filter((item) => item.celltype === group).map((item) => item.x)
    //
    //              trace.y.push(...values) // 将该 cluster 的所有值添加到 y 数组
    //              trace.x.push(...Array(values.length).fill(group)) // 将对应数量的 cluster 添加到 x 数组
    //
    //              return trace
    //            })
    //
    //            //  console.log(newboxPlotData)
    //            this.boxPlotData = newboxPlotData
    //            this.isLoading = false
    //            this.loading = true
    //            axios
    //              .post('http://scisd.databio1.com/api/query_diffgene', {
    //                subject: this.selectedGSE,
    //                cluster: this.celltype.replace(/\s+/g, '_')
    //              })
    //              .then((response1) => {
    //                //  console.log(response1.data)
    //                this.tabledata = response1.data
    //                this.imgurlx =
    //                  'http://scisd.databio1.com/api/cellmarker_enrichment_BP/' +
    //                  this.selectedGSE +
    //                  '/' +
    //                  this.celltype.replace(/\s+/g, '_')
    //                axios
    //                  .get(
    //                    'http://scisd.databio1.com/api/cellmarker_enrichment_KEGG_table/' +
    //                      this.selectedGSE +
    //                      '/' +
    //                      this.celltype.replace(/\s+/g, '_')
    //                  )
    //                  .then((response2) => {
    //                    console.log(JSON.parse(response2.data))
    //                    this.tabledata1 = JSON.parse(response2.data)
    //                    //get genelist
    //                    axios
    //                      .get('http://scisd.databio1.com/api/gene_list/' + this.selectedGSE)
    //                      .then((responseb) => {
    //                        const x = JSON.parse(responseb.data)
    //                        const geneArray1 = x.map((item) => Object.values(item)[0])
    //                        this.genes = geneArray1
    //
    //                        axios
    //                          .get('http://scisd.databio1.com/api/vs_table/' + this.selectedGSE)
    //                          .then((responseb) => {
    //                            const x = JSON.parse(responseb.data)
    //                            const geneArray1 = x.map((item) => Object.values(item)[1])
    //                            //  this.genes = geneArray1
    //                            this.multigroups = geneArray1
    //                            console.log(x)
    //                          })
    //                          .catch((error) => {
    //                            console.error(error)
    //                          })
    //                      })
    //                      .catch((error) => {
    //                        console.error(error)
    //                      })
    //                  })
    //                  .catch((error) => {
    //                    console.error(error)
    //                  })
    //              })
    //              .catch((error) => {
    //                console.error(error)
    //              })
    //          })
    //          .catch((error) => {
    //            console.error(error)
    //          })
    //      }
    //    }
  }
}
</script>
