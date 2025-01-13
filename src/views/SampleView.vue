<template>
  <div class="container mx-auto">
    <div class="card">
      <DataTable :value="tabledata" tableStyle="min-width: 50rem" stripedRows>
        <Column field="Subject" header="Subject">
          <template #body="slotProps">
            <router-link
              v-if="slotProps.data.Subject"
              :to="{ name: 'results', params: { id: slotProps.data.Subject } }"
              class="underline text-blue-500 hover:text-blue-700"
            >
              {{ slotProps.data.Subject }}
            </router-link>
          </template>
        </Column>
        <Column field="PMID" header="PMID"></Column>
        <Column field="Sample_count" header="Sample Counts"></Column>

        <Column field="Sample_disease_lesional" header="Lesional Samples"></Column>
        <Column field="Sample_disease_nonlesional" header="Nonlesional Samples"></Column>
        <Column field="Sample_healthy_control" header="Healthy Controls"></Column>
        <Column field="Sample_source" header="Sample Source"></Column>

        <Column field="disease_type" header="Disease Type"></Column>
        <Column field="n_cells" header="Cell Counts"></Column>
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
    const response = await axios.post('http://scisd.databio1.com/api/summary_database', {
      subject: 'all'
    })
    tabledata.value = response.data
    //  console.log(response.data)
  } catch (error) {
    console.error(error)
  }
})
</script>
