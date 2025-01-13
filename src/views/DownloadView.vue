<template>
  <div class="container mx-auto p-4">
    <div class="card">
      <DataTable :value="tabledata" tableStyle="min-width: 50rem" stripedRows>
        <Column field="subject" header="Subject"> </Column>
        <Column field="descript" header="Cell numbers">
          <template #body="slotProps">
            {{ slotProps.data.descript[0].n_cells }}
          </template>
        </Column>
        <Column field="gene_count" header="Gene numbers"> </Column>
        <Column field="file_size" header="File size"> </Column>
        <Column field="descript" header="Disease type">
          <template #body="slotProps">
            {{ slotProps.data.descript[0].disease_type }}
          </template>
        </Column>
        <Column field="descript" header="PMID">
          <template #body="slotProps">
            <a
              :href="`https://pubmed.ncbi.nlm.nih.gov/${slotProps.data.descript[0].PMID}`"
              class="underline text-blue-500 hover:text-blue-700"
            >
              {{ slotProps.data.descript[0].PMID }}
            </a>
          </template>
        </Column>
        <Column field="download" header="Download">
          <template #body="slotProps">
            <a
              :href="`http://scisd.databio1.com/api${slotProps.data.download_url}`"
              class="underline text-blue-500 hover:text-blue-700"
            >
              SeuratObject.rds
            </a>
          </template>
        </Column>
      </DataTable>
    </div>
  </div>
</template>

<script setup>
import { ref, onMounted } from 'vue'
import DataTable from 'primevue/datatable'
import Column from 'primevue/column'
import axios from 'axios'

const tabledata = ref([])

onMounted(async () => {
  try {
    const response = await axios.get('http://scisd.databio1.com/api/list_files')
    tabledata.value = response.data
    //console.log(response.data)
  } catch (error) {
    console.error(error)
  }
})
</script>
